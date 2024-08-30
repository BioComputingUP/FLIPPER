#!/usr/bin/env python
import math
import os
import logging.config
import configparser
import argparse
import gzip
import tempfile
import shutil
import copy
import sys
import json
import re

import numpy as np
import pandas as pd
from Bio.SVDSuperimposer import SVDSuperimposer

from Bio.PDB import ShrakeRupley
from Bio.PDB.MMCIFParser import MMCIF2Dict, FastMMCIFParser
from Bio.SeqUtils import IUPACData

from Bio.PDB.Polypeptide import is_aa, PPBuilder

import flipper

# silence scikit unpickle warning
from warnings import simplefilter
simplefilter(action='ignore', category=UserWarning)


# Secondary structure ranges in
RAMA_SS_RANGES = [(-180, -180, 80, 60, 'E', 'blue'),
                  (-180, 50, 80, 130, 'E', 'blue'),
                  (-100, -180, 100, 60, ' ', 'green'),  # poly-proline
                  (-100, 50, 100, 130, ' ', 'green'),  # poly-proline
                  (-180, -120, 180, 170, 'H', 'red'),
                  (0, -180, 180, 360, ' ', 'yellow')]  # loop

# Sander: Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
RSA_SCALE = {
    "ALA": 106.0,
    "ARG": 248.0,
    "ASN": 157.0,
    "ASP": 163.0,
    "CYS": 135.0,
    "GLN": 198.0,
    "GLU": 194.0,
    "GLY": 84.0,
    "HIS": 184.0,
    "ILE": 169.0,
    "LEU": 164.0,
    "LYS": 205.0,
    "MET": 188.0,
    "PHE": 197.0,
    "PRO": 136.0,
    "SER": 130.0,
    "THR": 142.0,
    "TRP": 227.0,
    "TYR": 222.0,
    "VAL": 142.0,
}


# resolution : [ Wilson_B, Average_B ]
WILSON_B = [
 [0.00, 10.0, 13.11],
 [1.00, 11.0, 16.44],
 [1.25, 14.0, 19.14],
 [1.50, 18.0, 21.76],
 [1.75, 23.0, 26.82],
 [2.25, 36.0, 39.42],
 [2.50, 44.0, 44.73],
 [2.75, 54.0, 51.94],
 [3.00, 66.0, 60.76],
 [3.25, 82.0, 78.70],
 [3.50, 93.0, 88.84],
 [3.75, 112.0, 102.29],
 [4.00, 135.0, 121.349],
 [4.25, 162.0, 143.960],
 [4.50, 194.0, 170.784],
 [4.75, 233.0, 202.606],
 [5.00, 280.0, 240.357],
 [5.25, 336.0, 285.142],
 [5.50, 404.0, 338.272],
 [5.75, 485.0, 401.301],
 [6.00, 550.0, 550.00]
]


def get_structure(pdb_file_gz):

    # Extract the PDB name. Useful for the log
    pdb_id = os.path.basename(pdb_file_gz).split(".")[0]  # Assume file base name is pdb2zpm.ent.gz

    # Uncompress the pdb file. The extension is important for FLIPPER
    pdb_file = "{}/{}.cif".format(tmp_dir, pdb_id)
    with open(pdb_file, 'w') as fout:
        with gzip.open(pdb_file_gz, 'rb') as f:
            for line in f:
                fout.write(line.decode())

    # Parse the original PDB
    try:
        # structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)
        structure = FastMMCIFParser(auth_chains=False, auth_residues=True, QUIET=True).get_structure(pdb_id, pdb_file)
    except Exception as err:
        logging.error("{} failed parsing cif file {}. ERROR: {} ".format(pdb_id, pdb_file, err))
        sys.exit(0)
    else:
        try:
            mmcif_dict = MMCIF2Dict(pdb_file)
        except Exception as err:
            logging.error("{} failed parsing cif dict {}. ERROR: {} ".format(pdb_id, pdb_file, err))
        else:
            res = mmcif_dict.get('_refine.ls_d_res_high') if mmcif_dict.get('_refine.ls_d_res_high') != ['.'] else mmcif_dict.get('_em_3d_reconstruction.res')  # 2.5
            if res is not None:
                res = float(res[0])

            met = mmcif_dict.get('_exptl.method')  # ["X-RAY DIFFRACTION"]
            if met is not None:
                met = met[0]

        logging.debug(
            "{} start, resolution {}, method {}, number models {}, number chains {}".format(pdb_id, res,
                                                                                                         met,
                                                                                                         len(structure),
                                                                                                         len(structure[0])))
        return pdb_file, pdb_id, structure, res, met


def get_secondary_structure(structure):
    """
    PPBuilder method doesn't work with CA only chains because can't calculate PHI and PSI
    """
    data = []  # [model_id, chain_id, het_code, pos, icode, index, aa, ss, phi, psi, phi_deg, psi_deg]
    ppb = PPBuilder()  # PolyPeptideBuilder
    for j, model in enumerate(structure):
        for chain in model:
            # Check there are residues in the chain
            if list(filter(is_aa, model[chain.id].get_residues())):
                c = 0
                for pp in ppb.build_peptides(model[chain.id]):
                    phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
                    for i, residue in enumerate(pp):
                        phi, psi = phi_psi[i]
                        phi_ = math.degrees(phi) if phi else None
                        psi_ = math.degrees(psi) if psi else None
                        # resname = IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())
                        resname = residue.get_resname()

                        ss_c_ = None
                        # Identify the SS class according to PSI, PHI and Ramachandran regions
                        if phi is not None and psi is not None:
                            for x, y, width, height, ss_c, color in RAMA_SS_RANGES:
                                if x <= phi_ < x + width and y <= psi_ < y + height:
                                    ss_c_ = ss_c
                                    break

                        data.append([model.serial_num, chain.id, *residue.get_full_id()[3], c, resname, ss_c_, phi, psi, phi_, psi_])
                        c += 1
            else:
                logging.debug("Zero valid AAs in chain {}".format(chain.id))

    return pd.DataFrame(data, columns=["model_id", "chain_id", "het", "pos", "ins", "index", "aa", "ss", "phi", "psi", "phi_deg", "psi_deg"])


def get_asa(structure, n_points=100):

    data = []  # [chain_id, het_code, pos, icode, index, aa, asa, rsa]
    data_isolation = []  # [chain_id, het_code, pos, icode, index, aa, asa, rsa]

    model_id = structure[0].serial_num

    # Calculate SASA for the entire complex
    sr = ShrakeRupley(n_points=n_points)
    sr.compute(structure[0], level='R')

    for chain in structure[0]:
        for i, residue in enumerate(structure[0][chain.id].get_residues()):
            if is_aa(residue):
                rsa_norm = RSA_SCALE.get(residue.get_resname())
                # resname = IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())
                resname = residue.get_resname()
                if rsa_norm:
                    data.append([model_id, chain.id, *residue.get_full_id()[3], i, resname, residue.sasa, residue.sasa / rsa_norm if rsa_norm else None])

    # If multichain calculates delta-rsa between the complexed chain and the chain in isolation
    if len(structure[0]) > 1:
        for chain in structure[0]:

            # Create a copy of the chain in isolation
            struct_chain = copy.deepcopy(structure[0][chain.id])
            sr.compute(struct_chain, level='R')

            for i, residue in enumerate(struct_chain.get_residues()):
                if is_aa(residue):

                    rsa_norm = RSA_SCALE.get(residue.get_resname())
                    # resname = IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())
                    resname = residue.get_resname()
                    if rsa_norm:
                        data_isolation.append([model_id, chain.id, *residue.get_full_id()[3], i, resname, residue.sasa, residue.sasa / rsa_norm if rsa_norm else None])

    columns = ["model_id", "chain_id", "het", "pos", "ins", "index", "aa", "sasa", "rsa"]
    return pd.merge(pd.DataFrame(data, columns=columns), pd.DataFrame(data_isolation, columns=columns), how='left', on=["model_id", "chain_id", "het", "pos", "ins", "index", "aa"], suffixes=["", "_isolation"])


def flipper_dicts(df, pdb_id, column="sasa"):
    data = {}
    for i, df_g in df.groupby("model_id"):
        for index, row in df_g.iterrows():
            if ~row.isnull().any():
                # "{}_{}_{}_{}{}".format(self.pdb_id, self.model_id, self.chain_id, self.pdb_index, self.pdb_insertion_code)
                residue_id = "{}_0_{}_{}{}".format(pdb_id, row["chain_id"], row["pos"], row["ins"] if row["ins"] != " " else " ")  # always set model_id to 0
                data[residue_id] = (row["ss"], row[column])
        if i > 0:
            break
    return data


def get_bfactor(structure, resolution, config_params):

    bfactor = []
    for chain in structure[0]:
        c = 0
        for atom in structure[0][chain.id].get_atoms():
            if atom.get_name() == "CA":
                bfactor.append([chain.id, *atom.get_parent().get_full_id()[3], atom.get_parent().get_resname(), c, atom.get_bfactor()])
                c += 1

    if bfactor:
        df_bfactor = pd.DataFrame(bfactor, columns=["chain_id", "het", "pos", "ins", "aa", "index", "bfactor"])

        # Find the closest resolution in the Wilson table
        wb = np.array(WILSON_B)
        index = int(np.abs(wb[:, 0] - resolution).argmin())
        wilson_b, b_average = wb[index][0], wb[index][1]

        # Set disorder
        th = wilson_b * config_params.getfloat("wilson_b_factor")
        logging.debug("bfactor threshold: {}".format(th))

        df_bfactor["bfactor_normalized"] = df_bfactor["bfactor"] / (4.0 * th)
        df_bfactor.loc[df_bfactor["bfactor"] > (4.0 * th), "bfactor_normalized"] = 1.0
    else:
        return

    return df_bfactor


def get_flipper(pdb_id, pdb_file, sasa_dict, sasa_dict_chains, config_flipper):
    """
    Features (for each residue for each chain):
    0 inter contacts (window average)
    1 intra long range contacts (window average)
    2 helix
    3 beta
    4 non-ss
    5 rsa
    6 delta-rsa
    7 inter contacts
    8 long range contacts
    9 distance_3D (linearity)
    10 length cutoff
    """

    data = []
    columns = "chain_id", "pos", "ins", "lip", "lip_status", "inter_contacts", "intra_long_contacts", "helix", "beta", "coil", "sasa", "delta_sasa", "linearity", "length_cutoff", "inter_contacts_chains"

    predictions, raw_pred, struct, neighbors, features = flipper.run_flipper(pdb_id, pdb_file, sasa_dict, sasa_dict_chains, config_flipper)

    interactions = {}
    for residue_id in neighbors.nn:
        chain_id = residue_id.split("_")[2]
        interactions.setdefault(chain_id, []).append(list(set([res.chain_id for res in neighbors.nn[residue_id][0]])))

    for chain_id in predictions.keys():
        for i, (p, pr, feat, inter_chains) in enumerate(zip(predictions[chain_id], raw_pred[chain_id], features[chain_id], interactions[chain_id])):
            data.append([
                chain_id,
                struct.chains[chain_id].residues[i].pdb_index,
                struct.chains[chain_id].residues[i].pdb_insertion_code.strip(),
                round(float(pr), 3),
                str(p),
                round(feat[0], 3),
                round(feat[1], 3),
                round(feat[2], 3),
                round(feat[3], 3),
                round(feat[4], 3),
                round(feat[5], 3),
                round(feat[6], 3),
                round(feat[9], 3),
                round(feat[10], 3),
                inter_chains
            ])

    return pd.DataFrame(data, columns=columns)


def superimpose(segment_1, segment_2):
    """
    return the new coordinates of segment_2
    upon superposition with segment_1
    """
    sup = SVDSuperimposer()
    sup.set(segment_1, segment_2)
    sup.run()
    rot, tran = sup.get_rotran()
    return np.dot(segment_2, rot) + tran  # the updated coordinate is returned


def get_mobi_distance(structure, config_params):
    residues = []  # [[chain_id, het, pos, ins, resname, index], ...]
    coords = {}  # {chain_id: []}
    for i, model in enumerate(structure):
        for chain in model:
            c = 0
            atoms = []
            for atom in model[chain.id].get_atoms():
                if atom.get_name() == "CA":
                    atom_coords = atom.get_coord()
                    if atom_coords is not None:
                        atoms.append(atom_coords)
                    if i == 0:
                        residues.append([chain.id, *atom.get_parent().get_full_id()[3], atom.get_parent().get_resname(), c])
                        c += 1

            # Check coordinates are provided for all atoms (sse residue 52 of model 11 of 1lpv
            if atoms and len(atoms) == len(residues):
                coords.setdefault(chain.id, []).append(atoms)  # one list of atoms for each model
            else:
                logging.debug("No CA atoms in chain {}".format(chain.id))
    df_data = pd.DataFrame(residues, columns=["chain_id", "het", "pos", "ins", "aa", "index"])

    df_list = []  # list of chain dataframes
    for chain_id in coords:
        if len(coords[chain_id]) > 1:
            # Align models
            df_coords = pd.DataFrame(coords[chain_id])
            coords_list = []
            coords_list.append(np.asarray(df_coords.iloc[0, :].values.tolist()))
            for i in range(1, df_coords.shape[0]):
                # Calculate new coordinates
                coords_list.append(superimpose(coords_list[0], np.asarray(df_coords.iloc[i].values.tolist())))

            # Calculate RMSD at the residue level (all Vs first model)
            res_distance = []
            # Swap axes
            coords_array = np.array(coords_list)
            coords_array = np.swapaxes(coords_array, 0, 1)
            for i in range(coords_array.shape[0]):
                distance = []
                for j in range(1, coords_array.shape[1]):
                    distance.append(np.sqrt(np.sum(coords_array[i][0] - coords_array[i][j]) ** 2))
                res_distance.append(distance)
            res_distance = np.array(res_distance)

            # Calculate Mobi scaled distance
            scaled_distance = []
            for distance in res_distance:
                # tmp_dist = 1.0 / (1.0 + pow((ca - ca_list_ref[j]) / config_params.getfloat('mobi_d_0'), 2.0))
                scaled_distance.append(1.0 / (1.0 + pow(distance / config_params.getfloat('mobi_d_0'), 2)))
            scaled_distance = np.array(scaled_distance)

            scaled_distance_mean = scaled_distance.mean(axis=1)
            scaled_distance_std = scaled_distance.std(axis=1)

            # Combine mena and std into a dataframe
            df_scaled_distance = pd.DataFrame(np.array([scaled_distance_mean, scaled_distance_std]).T, columns=["scaled_distance_mean", "scaled_distance_std"]).reset_index()
            df_scaled_distance['chain_id'] = chain_id

            df_list.append(pd.DataFrame.merge(df_data[df_data["chain_id"]==chain_id], df_scaled_distance, on=["chain_id", "index"], how="left"))
        else:
            logging.debug("Only one or zero models in chain {}".format(chain_id))
    if df_list:
        return pd.concat(df_list)  # Concatenate dataframes of each chain
    else:
        return None


def get_mobi_state(mobi_state, mobi_filter):
    # Filter for patterns
    mobi_patterns = [('1011', '1111'), ('1101', '1111'), ('10011', '11111'), ('11001', '11111'),('01010', '00000'), ('00100', '00000'), ('001100', '000000')]
    mobile_str = "".join(map(str, mobi_state.astype(int)))
    mobile_filter = "".join(map(str, mobi_filter.astype(int)))

    for ori, rep in mobi_patterns:
        mobile_str = mobile_str.replace(ori, rep)

    # Further filtering for patterns
    for pattern in ['110', '011']:
        for m in re.finditer(pattern, mobile_str):
            pos = m.start()
            if pattern == '110':
                if mobile_filter[pos + 2]:
                    mobile_str = mobile_str[0:pos] + '111' + mobile_str[pos + 3:]
            else:
                if mobile_filter[pos]:
                    mobile_str = mobile_str[0:pos] + '111' + mobile_str[pos + 3:]
    return mobile_str


def get_mobi(df_mobi_distance, df_ss, config_params):

    # Calculate secondary structure variation across models
    df_ = df_ss.groupby(["chain_id", "het", "pos", "ins", "aa", "index"]).agg(ss_same=("ss", lambda x: np.all(x.to_numpy() == x.iloc[0])),
                                                                              ss_first=("ss", "first"),
                                                                              phi_deg=("phi_deg", "mean"),
                                                                              psi_deg=("psi_deg","mean"))
    df_ = pd.merge(df_, df_mobi_distance, how="left", on=["chain_id", "het", "pos", "ins", "aa", "index"])

    # Mobile score is based on distance mean
    df_["mobi_score"] = 1 - df_["scaled_distance_mean"]

    # Secondary structure must be same and different from disorder " ".
    df_["ss_"] = ~df_["ss_same"] | (df_["ss_first"] == " ")

    # Apply cutoffs
    df_["phi_deg_"] = df_["phi_deg"] > config_params.getfloat('mobi_phi')
    df_["psi_deg_"] = df_["psi_deg"] > config_params.getfloat('mobi_psi')
    df_["scaled_distance_std_"] = df_["scaled_distance_std"] > config_params.getfloat('mobi_d_std')

    # Combine scaled distance mean with SS judgement.
    # If secondary structure vary (ss_ == True), only then apply scaled distance cutoff
    mobi_filter = df_["psi_deg_"] & df_["phi_deg_"] & df_["scaled_distance_std_"]
    mobi_state = (df_["scaled_distance_mean"] < config_params.getfloat('mobi_d_mean')) & df_["ss_"]
    df_["mobi_state"] = list(get_mobi_state(mobi_state, mobi_filter))

    columns = ["chain_id", "het", "pos", "ins", "aa", "index", "mobi_score", "mobi_state"]

    return df_[columns]


def write_doc(pdb_id, df_out, out_file):
    # {"pdb_id": "1jsu", "chain_id": "A", "residue_id": "13", "dssp": "-", "rsa": 0.54,
    # "bfactor": 70.19, "bfactor_normalized": 1.0, "lip": 0.712, "lip_status": "0", "inter_contacts": 2.167,
    # "intra_long_contacts": 0.167, "helix": 0.0, "beta": 0.333, "coil": 0.667, "delta_rsa": 0.224,
    # "linearity": 0.898, "length_cutoff": 1.0, "inter_contacts_chains": ["C"]}

    df_out = df_out.rename(columns={"chain_id": "label_asym_id", "pos": "pdb_residue_id", "inter_contacts_chains": "inter_contacts_label_asym_id"})

    df_out["pdb_id"] = pdb_id
    if "ss" in df_out:
        df_out["dssp"] = df_out["ss"]  # Fake DSSP. Based on PHI and PSI
        df_out.loc[df_out["dssp"]==" ", "dssp"] = "-"

    df_out["pdb_residue_id"] = df_out["pdb_residue_id"].map(str) + df_out["ins"]

    columns = ["pdb_id", "label_asym_id", "pdb_residue_id", "rsa", "bfactor", "bfactor_normalized", "lip", "lip_status",
               "inter_contacts", "intra_long_contacts", "dssp", "helix", "beta", "coil", "delta_rsa", "linearity",
               "length_cutoff", "inter_contacts_label_asym_id", "mobi_state", "mobi_score"]

    # Write document
    with gzip.open(out_file, "wb") as fout:
        for i, row in df_out.loc[:, df_out.columns.isin(columns)].iterrows():
            fout.write((row.dropna().to_json() + "\n").encode())
    return


def parse_args():
    parser = argparse.ArgumentParser(prog='biodb_disorder.py',
                                     description="Generate a JSON document with disorder annotation from a mmCIF file",
                                     epilog="Example: python biodb_disorder.py 2zps.cif.gz 2zps_disorder.mjson.gz")
    parser.add_argument('input_file', type=str, help="A gzip mmCIF file")
    parser.add_argument("output_file", type=str, help="Output JSON file. With \".gz\" extention")
    parser.add_argument("-ll", "--log_level", type=str,
                        choices=["notset", "debug", "info", "warning", "error", "critical"], default="info",
                        help="The log level")

    return parser.parse_known_args()


if __name__ == "__main__":
    # parse command line arguments
    args, unknown = parse_args()

    # Set logger
    logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s',
                        level=logging.getLevelName(args.log_level.upper()), stream=sys.stdout)

    # Setting configuration from file
    config = configparser.ConfigParser()
    config_file = "{}/config.ini".format(os.path.dirname(os.path.realpath(__file__)))
    config.read(config_file)
    logging.debug("Config file: {}".format(config_file))
    config = config["DEFAULT"]

    config_file_flipper = os.path.dirname(os.path.realpath(__file__)) + "/config_flipper.json"
    logging.debug("Config file flipper: {}".format(config_file_flipper))

    with open(config_file_flipper) as f:
        config_flipper = json.load(f)

    # Set temporary file names
    tmp_dir = tempfile.mkdtemp(prefix="disorder_")

    pdb_file, pdb_id, struct, res, met = get_structure(args.input_file)
    logging.info("{} pdb_id parsed, now processing".format(pdb_id))

    df_out = pd.DataFrame()

    # Bfactor
    if res:
        df_bfactor = get_bfactor(struct, res, config)
        df_out = pd.concat([df_out, df_bfactor])
        logging.debug("{} bfactor calculated".format(pdb_id))

    # Secondary structure is provided for all models
    df_ss = get_secondary_structure(struct)

    if not df_ss.empty:
        logging.debug("{} SS calculated".format(pdb_id))

        if not df_out.empty:
            df_out = pd.merge(df_out, df_ss.loc[df_ss["model_id"]==df_ss.iloc[0, 0]], on=["chain_id", "het", "pos", "ins", "aa", "index"])
        else:
            df_out = df_ss.loc[df_ss["model_id"]==df_ss.iloc[0,0]]

        # MOBI
        if len(struct) > 1:
            df_mobi_distance = get_mobi_distance(struct, config)
            if df_mobi_distance is not None and not df_mobi_distance.empty:
                df_mobi = get_mobi(df_mobi_distance, df_ss, config)
                if not df_out.empty:
                    df_out = pd.merge(df_out, df_mobi, on=["chain_id", "het", "pos", "ins", "aa", "index"])
                else:
                    df_out = df_mobi
                logging.debug("{} mobi calculated".format(pdb_id))

            else:
                logging.debug("Scaled distance not available can't calculate MOBI: {}".format(pdb_id))

        # LIPs prediction (FLIPPER)
        if len(struct[0]) > 1:
            # RSA is provided for all chains, both in complex and isolation but only for the first model
            df_sasa = get_asa(struct, n_points=(100 if len(struct[0]) < 10 else 10))
            if not df_sasa.empty:
                logging.debug("{} sasa calculated".format(pdb_id))

                # Prepare FLIPPER input (ASA and SS)
                df_ = pd.merge(df_ss, df_sasa, how='left', on=["model_id", "chain_id", "het", "pos", "ins", "index", "aa"])
                flipper_dict = flipper_dicts(df_, pdb_id, "rsa")
                flipper_dict_chains = flipper_dicts(df_, pdb_id, "rsa_isolation")
                # flipper_dict = flipper_dicts(df_, pdb_id, "sasa")
                # flipper_dict_chains = flipper_dicts(df_, pdb_id, "sasa_isolation")

                df_flipper = get_flipper(pdb_id, pdb_file, flipper_dict, flipper_dict_chains, config_flipper)
                df_flipper['ins'].replace("", " ", inplace=True)  # replace empty insertion codes
                if not df_out.empty:
                    df_out = pd.merge(df_out, df_flipper, on=["chain_id", "pos", "ins"])
                else:
                    df_out = df_flipper
                logging.debug("{} flipper calculated".format(pdb_id))

            else:
                logging.debug("SASA not available can't calculate FLIPPER: {}".format(pdb_id))
    else:
        logging.debug("SS not available can't calculate MOBI and FLIPPER: {}".format(pdb_id))

    if not df_out.empty:
        write_doc(pdb_id, df_out, args.output_file)
    else:
        logging.debug("Nothing to write on output: {}".format(pdb_id))

    # Clean temporary files when completing
    shutil.rmtree(tmp_dir, ignore_errors=True)

    # ls ../biodb_data/cif/ | while read line; do python disorder/biodb_disorder.py ../biodb_data/cif/$line ../biodb_data/disorder_$line; done
