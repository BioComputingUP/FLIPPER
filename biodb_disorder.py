#!/usr/bin/env python

import os
import logging.config
import configparser
import argparse
import gzip
import tempfile
import shutil
import signal
import sys
import json

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import DSSP
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.Polypeptide import is_aa

import mobi
import flipper


class ModelSelect(Select):
    """
    Select model for BioPython PDB save
    """

    def __init__(self, model_ids):
        # self.model_id = model_id
        self.model_ids = model_ids

    def accept_model(self, model):
        # print self.model_id, model.serial_num
        # return (self.model_id + 1) == model.serial_num
        return (model.serial_num) in self.model_ids


class ChainSelect(Select):
    """
    Select model for BioPython PDB save
    """

    def __init__(self, chain_ids):
        # self.model_id = model_id
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        # print self.model_id, model.serial_num
        # return (self.model_id + 1) == model.serial_num
        return chain.id in self.chain_ids


class NotAltLoc(Select):
    def accept_atom(self, atom):
        return atom.get_altloc() == ' '


def write_models(pdb_id, structure, config_params):
    # Write temporary PDB model files. One different file for each chain and model
    # Models can start from a number greater than 0

    chain_model = {}  # {chain: [(model_id, model_file, model_dssp), ...]}
    model_chain = {}  # {model_id: [(chain_id, model_file, model_dssp), ...]}

    io = PDBIO()

    # Split models
    model_files = [] # [(model_id, model_file), ...]
    for model in structure:
        i = model.serial_num
        model_file = "{}/model_{}_{}".format(tmp_dir, pdb_id, i)
        model_files.append((i, model_file))

        # Save model
        io.set_structure(structure)
        io.save(model_file, select=ModelSelect([i]))

    # Split chains. For each model save chain separately and get dssp
    for model_id, model_file in model_files:

        model_structure = PDBParser(QUIET=True).get_structure(model_id, model_file)

        # model_structure should have length 1
        for model in model_structure:
            # print(i, model_file, model_structure, len(model_structure))
            for chain in model:
                # Save chain models separately
                model_chain_file = "{}_{}".format(model_file, chain.id)

                io.set_structure(model_structure)
                io.save(model_chain_file, select=ChainSelect([chain.id]))

                # Calculate secondary structure for each chain separately (for Flipper)
                model_structure_chain = PDBParser(QUIET=True).get_structure("{}_{}_{}".format(pdb_id, model_id, chain.id), model_chain_file)
                _dssp = None
                try:
                    _dssp = DSSP(model_structure_chain[0], model_chain_file, dssp=config_params.get('dssp'))
                except Exception as e:
                    logging.warning("{} DSSP error {} {} {}".format(pdb_id, model_file, chain.id, e))

                chain_model.setdefault(chain.id, []).append((model_id, model_chain_file, _dssp))
                model_chain.setdefault(model_id, []).append((chain.id, model_chain_file, _dssp))

    return chain_model, model_chain


def get_disorder(pdb_file_gz, out_file, config_params, config_flipper):
    docs = {}

    # Extract the PDB name. Useful for the log
    pdb_id = os.path.basename(pdb_file_gz)[3:7]  # Assume file base name is pdb2zpm.ent.gz

    # Uncompress the pdb file
    pdb_file = "{}/{}.pdb".format(tmp_dir, pdb_id)
    with open(pdb_file, 'w') as fout:
        with gzip.open(pdb_file_gz, 'rb') as f:
            for line in f:
                fout.write(line.decode())

    # Parse the original PDB
    try:
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)
    except Exception:
        logging.error("{} failed parsing pdb file {} ".format(pdb_id, pdb_file))
        shutil.rmtree(tmp_dir)
    else:

        resolution = float(structure.header['resolution']) if structure.header['resolution'] else None
        method = structure.header['structure_method'].upper()

        logging.info(
            "{} start, resolution {}, method {}, number models {}, number chains {}".format(pdb_id, resolution,
                                                                                                         method,
                                                                                                         len(structure),
                                                                                                         len(structure[0])))

        # Get rsa and dssp
        # chain_model --> {chain: [(model_id, model_file, model_dssp), ...]}
        # model_chain --> {model_id: [(chain_id, model_file, model_dssp), ...]}
        chain_model, model_chain = write_models(pdb_id, structure, config_params)

        # Assign secondary structure and rsa to the output document
        if model_chain:
            logging.info("{} dssp success".format(pdb_id))

            for model_id in model_chain:
                for chain_id, model_file, model_dssp in model_chain[model_id]:
                    if model_dssp:
                        structure = PDBParser(QUIET=True).get_structure(pdb_id, model_file)
                        for chain in structure[0]:  # do not know why but it start from one even if different
                            for residue in chain:
                                if is_aa(residue):
                                    key = (chain.id, residue.id)
                                    residue_id_pdb = "{}{}".format(residue.id[1], residue.id[2] if residue.id[2] != " " else "")
                                    if key in model_dssp:
                                        _, _, ss, rsa, _, _, _, _, _, _, _, _, _, _ = model_dssp[key]
                                        # print(key, model_dssp[key])
                                        docs.setdefault(chain.id, {}).setdefault(residue_id_pdb, {})
                                        docs[chain.id][residue_id_pdb]["dssp"] = ss
                                        if rsa != "NA":
                                            docs[chain.id][residue_id_pdb]["rsa"] = round(rsa, 3)
                break

        # get bfactor
        # chain_model --> {chain: [(model_id, model_file, model_dssp), ...]}
        # model_chain --> {model_id: [(chain_id, model_file, model_dssp), ...]}
        if resolution and model_chain:
            logging.info("{} running bfactor".format(pdb_id))

            for model_id in model_chain:
                for chain_id, model_file, model_dssp in model_chain[model_id]:
                    structure = PDBParser(QUIET=True).get_structure(pdb_id, model_file)

                    residues, sequence, ca_list = mobi._get_ca_list(structure[0][chain_id])

                    if ca_list:
                        bfactor, bfactor_normalized = mobi.get_bfactor(ca_list, resolution,
                                                                       config_params.getfloat("wilson_b_factor"))
                        if bfactor_normalized is not None:
                            for residue, b, bn in zip(residues, bfactor, bfactor_normalized):
                                residue_id_pdb = "{}{}".format(*residue)
                                docs.setdefault(chain_id, {}).setdefault(residue_id_pdb, {})
                                docs[chain_id][residue_id_pdb]["bfactor"] = round(b, 3)
                                docs[chain_id][residue_id_pdb]["bfactor_normalized"] = round(bn, 3)
                        else:
                            logging.warning("{} {} {} bfactor failed".format(pdb_id, chain.id, model_id))
                    else:
                        logging.info("{} {} {} no bfactor".format(pdb_id, chain_id, model_id))
                break

        # Filter models without dssp
        # chain_model --> {chain: [(model_id, model_file, model_dssp), ...]}
        # model_chain --> {model_id: [(chain_id, model_file, model_dssp), ...]}
        # get mobility
        if len(model_chain) > 1:
            logging.info("{} running mobi".format(pdb_id))
            for chain_id in chain_model:
                ca_list = None
                residues = None
                for model_id, model_file, model_dssp in chain_model[chain_id]:
                    if model_dssp:
                        structure = PDBParser(QUIET=True).get_structure(pdb_id, model_file)
                        residues, sequence, ca_list = mobi._get_ca_list(structure[0][chain_id])
                        break

                if ca_list:
                    mobile_score, mobile_str = mobi.get_mobi(pdb_id, chain_id, chain_model[chain_id], config_params)
                    if mobile_score is not None:
                        for residue, m, ms in zip(residues, mobile_score, mobile_str):
                            residue_id_pdb = "{}{}".format(*residue)
                            docs.setdefault(chain_id, {}).setdefault(residue_id_pdb, {})
                            docs[chain_id][residue_id_pdb]["mobile"] = round(m, 3)
                            docs[chain_id][residue_id_pdb]["mobile_status"] = ms
                    else:
                        logging.warning("{} {} mobi failed".format(pdb_id, chain_id))
                else:
                    logging.info("{} {} no dssp for this chain".format(pdb_id, chain_id))

        # Flipper predictor
        # chain_model --> {chain: [(model_id, model_file, model_dssp), ...]}
        # model_chain --> {model_id: [(chain_id, model_file, model_dssp), ...]}
        if len(chain_model) > 1:
            logging.info("{} running flipper".format(pdb_id))

            structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)

            # flipper needs to calculate asa in bound state
            dssp_dict = None
            try:
                dssp_dict = DSSP(structure[0], pdb_file, dssp=config_params.get('dssp'))
            except Exception as e:
                logging.warning("{} DSSP error {}".format(pdb_id, e))

            if dssp_dict:
                for model_id in model_chain:
                    dssp_dict_chains = {}
                    for chain_id, model_file, model_dssp in model_chain[model_id]:
                        if model_dssp:
                            dssp_dict_chains.setdefault(chain_id, []).append(model_dssp)
                    if dssp_dict_chains:
                        predictions, raw_pred, struct, neighbors, features = flipper.run_flipper(pdb_id, pdb_file, dssp_dict, dssp_dict_chains, model_id, config_flipper)
                        """
                        Features (foe each residue for each chain):
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

                        interactions = {}
                        for residue_id in neighbors.nn:
                            # if neighbors.nn[residue_id][0]:
                            #     for res in neighbors.nn[residue_id][0]:
                            #         print(residue_id, res.chain_id)
                            chain_id = residue_id.split("_")[2]
                            interactions.setdefault(chain_id, []).append(list(set([res.chain_id for res in neighbors.nn[residue_id][0]])))

                        for chain_id in predictions.keys():
                            for i, (p, pr, feat, inter_chains) in enumerate(zip(predictions[chain_id], raw_pred[chain_id], features[chain_id], interactions[chain_id])):

                                residue_id_pdb = "{}{}".format(struct.chains[chain_id].residues[i].pdb_index,
                                                               struct.chains[chain_id].residues[i].pdb_insertion_code.strip())
                                docs.setdefault(chain_id, {}).setdefault(residue_id_pdb, {})
                                docs[chain_id][residue_id_pdb]["lip"] = round(float(pr), 3)
                                docs[chain_id][residue_id_pdb]["lip_status"] = str(p)
                                docs[chain_id][residue_id_pdb]["inter_contacts"] = round(feat[0], 3)
                                docs[chain_id][residue_id_pdb]["intra_long_contacts"] = round(feat[1], 3)
                                docs[chain_id][residue_id_pdb]["helix"] = round(feat[2], 3)
                                docs[chain_id][residue_id_pdb]["beta"] = round(feat[3], 3)
                                docs[chain_id][residue_id_pdb]["coil"] = round(feat[4], 3)
                                docs[chain_id][residue_id_pdb]["rsa"] = round(feat[5], 3)
                                docs[chain_id][residue_id_pdb]["delta_rsa"] = round(feat[6], 3)
                                docs[chain_id][residue_id_pdb]["linearity"] = round(feat[9], 3)
                                docs[chain_id][residue_id_pdb]["length_cutoff"] = round(feat[10], 3)

                                if inter_chains:
                                    # print(inter_chains)
                                    docs[chain_id][residue_id_pdb]["inter_contacts_chains"] = inter_chains

                    else:
                        logging.warning("{} {} flipper failed, no chain dssp on model".format(pdb_id, model_id))
                    break
            else:
                logging.warning("{} flipper failed, no dssp".format(pdb_id))

        # Write document
        if docs:
            with gzip.open(out_file, "wb") as fout:
                for chain_id in docs:
                    for residue_id_pdb in docs[chain_id]:
                        doc = {"pdb_id": pdb_id, "chain_id": chain_id, "residue_id": residue_id_pdb}
                        doc.update(docs[chain_id][residue_id_pdb])
                        fout.write((json.dumps(doc) + "\n").encode())
    return


def parse_args():
    parser = argparse.ArgumentParser(prog='biodb_disorder.py',
                                     description="Execute DSSP and TM-score to generate a JSON document with disorder annotation.",
                                     epilog="  Example: python3 biodb_disorder.py pdb2zps.ent.gz 2zps.mjson.gz")
    parser.add_argument('input_file', type=str, help="A gzip PDB file")
    parser.add_argument("output_file", type=str, help="Output JSON file. With \".gz\" extention")
    parser.add_argument("-ll", "--log_level", type=str,
                        choices=["notset", "debug", "info", "warning", "error", "critical"], default="info",
                        help="The log level")
    parser.add_argument("--production", action='store_true', default=False,
                        help="Activate the production section in the configuration file")

    parser.add_argument("--keep_tmp_file", action='store_true', default=False,
                        help="Don't remove temporary files")

    return parser.parse_args()


def signal_handler(sig, frame):
    logging.error('Received signal: {}'.format(sig))
    shutil.rmtree(tmp_dir)
    exit('Received signal: {}'.format(sig))  # Allow to exit with a signal different from 0


if __name__ == "__main__":
    # parse command line arguments
    args = parse_args()

    # Set logger
    logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s',
                        level=logging.getLevelName(args.log_level.upper()), stream=sys.stdout)

    # Setting configuration from file
    config = configparser.ConfigParser()
    config_file = "{}/config.ini".format(os.path.dirname(os.path.realpath(__file__)))
    config.read(config_file)
    logging.debug("Config file: {}".format(config_file))

    config_file_flipper = os.path.dirname(os.path.realpath(__file__)) + (
        "/config_flipper.json" if not args.production else "/config_flipper_production.json")
    logging.debug("Config file flipper: {}".format(config_file_flipper))
    with open(config_file_flipper) as f:
        config_flipper = json.load(f)

    # Set temporary file names
    tmp_dir = tempfile.mkdtemp(prefix="disorder_")

    # Signal listeners (to delete files even when interrupted/terminated)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    # Print pdb_id on stderr (use only for debugging)
    # sys.stderr.write(args.input_file + '\n')

    # Run the program
    get_disorder(args.input_file, args.output_file, config["production"] if args.production else config["DEFAULT"],
                 config_flipper)

    # Clean temporary files when completing
    if not args.keep_tmp_file:
        shutil.rmtree(tmp_dir)


    # Examples of problematic PDBs

    # NMR 2kkw
    # XRAY which insertion codes 1cu4
    # PDB with chain breaks

    # awk '{print($1)}' ../particular-pdbs.dat | while read line; do python3 biodb_disorder.py /db/pdb/"${line:1:2}"/pdb"$line".ent.gz test_"$line".gz; done
