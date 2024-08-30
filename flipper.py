#!/usr/bin/env python

import logging
from Bio.SeqUtils import IUPACData

from lips_finder import LipsFinder

sander = {
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
    "VAL": 142.0
}

wilke = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLN": 225.0,
    "GLU": 223.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


def dssp_rsa_to_asa(dssp, pdb_id, model_id):


    dssp_new = {}  # {residue_id: (ss, asa)}
    for k in dssp.keys():

        _, aa, ss, rsa, _, _, _, _, _, _, _, _, _, _ = dssp[k]

        # residue_id = "{}_{}_{}_{}{}".format(pdb_id, model_id, k[0], k[1][1], k[1][2] if k[1][2] != " " else " ")
        residue_id = "{}_0_{}_{}{}".format(pdb_id, k[0], k[1][1], k[1][2] if k[1][2] != " " else " ")  # always set model_id to 0
        rsa_new = 0.0
        try:
            # Convert to ASA (BioPython uses Sander by default) and then normalize by Wilke
            rsa_new = rsa * sander.get(IUPACData.protein_letters_1to3.get(aa, 'XXX').upper(), 0.0) / wilke.get(IUPACData.protein_letters_1to3.get(aa, 'XXX').upper(), 0.0)
        except Exception:
            logging.debug("{} {} {} {} DSSP rsa error".format(pdb_id, model_id, aa, rsa))
        # print(rsa_new, rsa, max_acc.get(IUPACData.protein_letters_1to3.get(aa, 'XXX').upper(), 0.0))
        dssp_new[residue_id] = (ss, rsa_new)

    return dssp_new


def run_flipper(pdb_id, pdb_file, dssp_dict, dssp_dict_chains, config_params):
    """
    Execute FLIPPER
    :param pdb_id: the PDB identifier, 4 letter code
    :param pdb_file: path to the PDB file
    :param dssp_dict:
    :param dssp_dict_chains:
    :param model_id:
    :param config_params:
    :return:
    """

    finder = LipsFinder(config_params['finder'], dssp_dict, dssp_dict_chains)
    finder.model_in()

    predictions, raw_pred, structure, neighbors, features = finder.predict(pdb_id,
                                                                 pdb_file,
                                                                 proba=config_params["proba"],
                                                                 blur=config_params["blur"],
                                                                 gap=config_params["gap_fill"],
                                                                 threshold=config_params["proba_threshold"],
                                                                 model_id=0)
    return predictions, raw_pred, structure, neighbors, features

