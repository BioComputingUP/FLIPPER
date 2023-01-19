#!/usr/bin/env python


import logging.config
import subprocess
import shlex
import copy
import re

import numpy as np
from Bio.SeqUtils import IUPACData
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


# resolution : [ Wilson_B, Average_B ]
WILSON_B = [
 [0.00, [10.0, 13.11]],
 [1.00, [11.0, 16.44]],
 [1.25, [14.0, 19.14]],
 [1.50, [18.0, 21.76]],
 [1.75, [23.0, 26.82]],
 [2.25, [36.0, 39.42]],
 [2.50, [44.0, 44.73]],
 [2.75, [54.0, 51.94]],
 [3.00, [66.0, 60.76]],
 [3.25, [82.0, 78.70]],
 [3.50, [93.0, 88.84]],
 [3.75, [112.0, 102.29]],
 [4.00, [135.0, 121.349]],
 [4.25, [162.0, 143.960]],
 [4.50, [194.0, 170.784]],
 [4.75, [233.0, 202.606]],
 [5.00, [280.0, 240.357]],
 [5.25, [336.0, 285.142]],
 [5.50, [404.0, 338.272]],
 [5.75, [485.0, 401.301]],
 [6.00, [550.0, 550.00]]
]


def _get_ca_list(chain):
    """

    :param chain: The structure chain object
    :return:
    """
    ca_list = []  # [<ca_atom_object or None>, ...]
    residues = []  # [(<index>, <insertion_code>, <3_letter_residue_name>_upper>), ...]
    sequence = ""

    for residue in chain:
        if is_aa(residue):
            _, _, chain_id, res_id = residue.get_full_id()
            try:
                residues.append((res_id[1], res_id[2].strip(), residue.get_resname().upper()))
                sequence += IUPACData.protein_letters_3to1_extended.get(residue.get_resname().capitalize(), 'X')
                ca_list.append(residue['CA'])
            except KeyError:
                logging.warning("Failed to find CA in residue {}".format(residue.get_full_id()))

    return residues, sequence, ca_list


def get_bfactor(ca_list, resolution, wilson_b_factor):

    bfactor = []
    for ca in ca_list:
        if ca is not None:
            bfactor.append(ca.get_bfactor())
        else:
            bfactor.append(None)

    if not all(x is None for x in bfactor):

        # Find the closest resolution in the Wilson table
        # resolution = np.abs(np.array(WILSON_B.keys()) - pdb_complex.resolution).min()
        # wilson_b, b_average = WILSON_B[resolution]

        wb = np.array(WILSON_B)
        index = int(np.abs(wb[:, 0] - resolution).argmin())
        wilson_b, b_average = wb[index, 1]
        # logging.debug("{} {} {} {}".format(wb, index, wilson_b, b_average))

        # Set disorder
        th = wilson_b * wilson_b_factor
        logging.debug("bfactor threshold: {}".format(th))

        # A bfactor normalized grater than 0.5 is disordered (i.e. greater than 2.0 * th)
        bfactor_normalized = []
        for l in bfactor:
            if l is None:
                bfactor_normalized.append(None)
            else:
                if l > 4.0 * th:
                    bfactor_normalized.append(1.0)
                else:
                    bfactor_normalized.append(l / (4.0 * th))
        bfactor_normalized = np.asarray(bfactor_normalized, dtype=float)

        return bfactor, bfactor_normalized

    return None, None


def get_mobi(pdb_id, chain_id, chain_models, config_params):
    """
    Implement Mobi
    mobi_scaled_distance = 1 / (1 + (atom_distance / d0)**2)
    """
    ca_list_ref = None
    ca_list = None
    phi = []
    psi = []
    ss = []
    distance_ca = []  # Higher numbers correspond to closer atoms !!!
    data = []  # just for printing

    ppb = PPBuilder(10.0)  # Build the chains allowing for "holes" of max 10 A

    # Each model file contains a single chain
    first_model = None
    for i, (model_id, model_file, model_dssp) in enumerate(chain_models):
        if model_dssp:
            # logging.debug(
            #     "{} {} Mobi processing: {} {} {} {}".format(pdb_id, chain_id, i, model_id, model_file, model_dssp))

            structure = PDBParser(QUIET=True).get_structure(pdb_id, model_file)

            # Calculate phi and psi
            phi_psi_list = []
            pp_list = ppb.build_peptides(structure[0])  # , aa_only=False)  # TODO check it does not want the chain instead
            for pp in pp_list:
                phi_psi_list += pp.get_phi_psi_list()

            if pp_list:
                logging.debug(
                    "{} {} Mobi processing: {} {} {} {}".format(pdb_id, chain_id, first_model, i, model_id, model_file, model_dssp))
                if first_model is None:
                    first_model = i
                    # Initialize chain arrays
                    residues, sequence, ca_list = _get_ca_list(structure[0][chain_id])
                    ca_list_ref = ca_list

                # Structural alignment (TM-align) all against first model
                if i > first_model:

                    command = '{} {} {} -o {}_tmscore'.format(config_params.get('tm_score'), model_file, chain_models[first_model][1], model_file)
                    logging.debug("{} alignment command: {}".format(pdb_id, command))
                    try:
                        subprocess.check_output(shlex.split(command))
                    except subprocess.CalledProcessError:
                        # Note TMalign crashes when PDBs have negative residue numbering
                        logging.error("{} {} TMscore error".format(pdb_id, chain_id, command))
                    else:
                        # Modify the output file to make it readable by BioPython (add last 2 columns)
                        alignment_file = "{}_tmscore.pdb".format(model_file)
                        try:
                            aligned_structure = PDBParser(QUIET=True).get_structure(pdb_id, alignment_file)
                        except Exception as e:
                            logging.error("{} failed to parse aligned PDB models".format(pdb_id, alignment_file))
                            # traceback.print_exc()  # Print to stderr
                            continue
                        else:
                            _, _, ca_list = _get_ca_list(aligned_structure[0][chain_id])

                # Break when ppb.build_peptides skip residues
                if len(ca_list) != len(ca_list_ref):
                    logging.warning("{} length error ca_list_ref {}, ca_list {} {} ".format(pdb_id, len(ca_list_ref), len(ca_list), model_file))
                else:
                    # Assign matrix distances, psi, phi, ss
                    for j, ca in enumerate(ca_list):

                        try:
                            ss.append(model_dssp[chain_id, ca.get_parent().get_full_id()[3]][2])
                        except Exception as e:
                            logging.debug("{} dssp missing key {}".format(pdb_id, e))
                            ss.append("-")

                        # Distances are calculated only from second model
                        if i > first_model:
                            tmp_dist = 1.0 / (1.0 + pow((ca - ca_list_ref[j]) / config_params.getfloat('mobi_d_0'), 2.0))
                            # print(j, ca - ca_list_ref[j], pow((ca - ca_list_ref[j]) / config_params.getfloat('mobi_d_0'), 2.0), tmp_dist, ca.get_vector(), ca_list_ref[j].get_vector())
                            distance_ca.append(tmp_dist)

                if len(ca_list_ref) != len(phi_psi_list):
                    logging.warning("{} length error ca_list_ref {}, phi_psi_list {} {} ".format(pdb_id, len(ca_list_ref), len(phi_psi_list), model_file))
                else:
                    # Assign matrix distances, psi, phi, ss
                    for j, ca in enumerate(ca_list):
                        phi.append(abs(phi_psi_list[j][0] * 180.0 / np.pi) if phi_psi_list[j][0] else 0)
                        psi.append(abs(phi_psi_list[j][1] * 180.0 / np.pi) if phi_psi_list[j][1] else 0)

            else:
                logging.warning("{} {} model not parsed {}".format(pdb_id, chain_id, model_file))

    if ca_list_ref is None:
        logging.warning("{} {} no parsed models".format(pdb_id, chain_id))
        return None, None

    # Calculate PSI and PHI deviation
    phi_state = None
    psi_state = None
    if phi:
        phi = np.asarray(phi).reshape((-1, len(ca_list_ref)))
        phi_state = np.std(phi, axis=0)  # STD by column
        data.append(("phi_std", phi_state))
        phi_state = phi_state > config_params.getfloat('mobi_phi')

    if psi:
        psi = np.asarray(psi).reshape((-1, len(ca_list_ref)))
        psi_state = np.std(psi, axis=0)  # STD by column
        data.append(("psi_std", psi_state))
        psi_state = psi_state > config_params.getfloat('mobi_psi')


    ss = np.asarray(ss, dtype="str").reshape((-1, len(ca_list_ref)))
    # Print the secondary structure for each position and each model
    # logging.debug("dssp_out {} {} {}".format(pdb_id, chain_id, ";".join(["{}{},{}".format(r[0], r[1], "".join(s)) for r, s in zip(residues, ss.T)])))

    distance_ca = np.asarray(distance_ca).reshape((-1, len(ca_list_ref)))

    logging.debug(distance_ca)
    if distance_ca.shape[0] < 2:
        logging.error("{} {} not enough models compared".format(pdb_id, chain_id))
        return None, None

    # Standard deviation of MOBI distances by column
    distance_std_state = np.std(distance_ca, axis=0)
    data.append(("distance_STD", distance_std_state))
    distance_std_state = distance_std_state > config_params.getfloat('mobi_d_std')
    data.append(("distance_STD_state", distance_std_state))

    # Average Scaled Distance. Represents closeness in reality
    mobile_score = np.mean(distance_ca, axis=0)
    data.append(("distance_mean", mobile_score))

    # Set mobile state
    mobile_str = copy.copy(mobile_score)
    mobile_str = mobile_str < config_params.getfloat('mobi_d_mean')
    data.append(("distance_mean_state", mobile_str))

    # Calculate disorder pattern from secondary structure
    # Transform the matrix in 1/0 by comparing with the first model
    ss_state = np.apply_along_axis(lambda y: [1 if aa == y[0] else 0 for aa in y], 0, ss)

    # If secondary structure is not constant, or is constant but ss is "S" or "-" , than disorder (True)
    ss_state = np.logical_or(~ss_state.all(axis=0), np.any(np.isin(ss, ['S', '-']), axis=0))
    data.append(("ss_state", ss_state))

    # Filter by SS assignment
    mobile_str[~ss_state] = False
    data.append(("distance_mean_state_ss", mobile_str))

    # Filter for patterns
    mobi_patterns = [('1011', '1111'), ('1101', '1111'), ('10011', '11111'), ('11001', '11111'),
                     ('01010', '00000'), ('00100', '00000'), ('001100', '000000')]
    mobile_str = "".join(map(str, mobile_str.astype(int)))
    for ori, rep in mobi_patterns:
        mobile_str = mobile_str.replace(ori, rep)
    data.append(("distance_mean_pattern", mobile_str))

    # Further filtering for patterns
    for pattern in ['110', '011']:
        for m in re.finditer(pattern, mobile_str):
            pos = m.start()
            if pattern == '110':
                if phi_state is not None and psi_state is not None and \
                        phi_state[pos+2] and psi_state[pos+2] and distance_std_state[pos+2] and psi_state[pos+1]:
                    mobile_str = mobile_str[0:pos] + '111' + mobile_str[pos+3:]
            else:
                if phi_state is not None and psi_state is not None and \
                        phi_state[pos] and psi_state[pos] and distance_std_state[pos] and psi_state[pos+1]:
                    mobile_str = mobile_str[0:pos] + '111' + mobile_str[pos+3:]
    data.append(("mobi_state", mobile_str))

    mobile_score = 1 - mobile_score  # Transform the distance into a score (high score high disorder)
    data.append(("mobi_score", mobile_score))

    logging.debug("ss_state: {}\nmobile score: {}\nmobile state: {}".format(ss_state, mobile_score, mobile_str))

    # for k, a in data:
    #     print(k, list(a))

    return mobile_score, mobile_str
