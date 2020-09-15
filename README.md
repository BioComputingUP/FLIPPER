FLIPPER, Fast Linear Interacting Peptides Predictor
==============================================
Damiano Piovesan, Paolo Bonato and Silvio C.E. Tosatto

Version 1.0

Introduction
------------
This package include the FLIPPER predictor and the MOBI software to extract LIPs and mobile residues from
PDB structures. FLIPPER is a random forest whereas MOBI simply compare the geometry of different models 
(when available) of the same PDB.

MOBI algorithm is described here:

    Mobi 2.0: An improved method to define intrinsic disorder, mobility and linear binding regions in protein structures
    Piovesan, D., Tosatto, S.C.E.
    (2018) Bioinformatics, 34 (1), pp. 122-123
    https://pubmed.ncbi.nlm.nih.gov/28968795/
    
    MOBI: A web server to define and visualize structural mobility in NMR protein ensembles
    Martin, A.J.M., Walsh, I., Tosatto, S.C.E.
    (2010) Bioinformatics, 26 (22), art. no. btq537, pp. 2916-2917
    https://pubmed.ncbi.nlm.nih.gov/20861031/

Requirements
------------
* Python 3
* TM-score compiled from the C++ version (https://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.cpp), 
this important as other sources generate a slightly different output
* DSSP 2.2.1 (ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.2.1.tgz). Again, the version matters as the output format change

Usage
-----
The following command prints the help page on screen:

    python3 biodb_disorder.py -h


By default, FLIPPER searches the ``ext_bin/`` folder (that includes TM-score and DSSP
executables) in its own directory. So, FLIPPER can be executed by just
providing the PDB input and the output file:

    python3 biodb_disorder.py pdb2zps.ent.gz 2zps.mjson.gz

Configuration
-------------
Check config.ini and config_flipper.json for correct paths to executables. 
Provide absolute paths. 

Output
------
Output is provided in multiline JSON format, e.g. one JSON document per line. 
Each document corresponds to one residue.

    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "25", "dssp": "-", "rsa": 0.539, "bfactor": 41.49, "bfactor_normalized": 0.24, "lip": 1.0, "lip_status": "1", "inter_contacts"
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "26", "dssp": "-", "rsa": 0.529, "bfactor": 27.85, "bfactor_normalized": 0.161, "lip": 0.994, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "27", "dssp": "T", "rsa": 0.529, "bfactor": 24.64, "bfactor_normalized": 0.143, "lip": 0.993, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "28", "dssp": "T", "rsa": 0.528, "bfactor": 17.94, "bfactor_normalized": 0.104, "lip": 1.0, "lip_status": "1", "inter_contacts
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "29", "dssp": "-", "rsa": 0.53, "bfactor": 21.7, "bfactor_normalized": 0.126, "lip": 0.99, "lip_status": "1", "inter_contacts"
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "30", "dssp": "-", "rsa": 0.53, "bfactor": 23.33, "bfactor_normalized": 0.135, "lip": 1.0, "lip_status": "1", "inter_contacts"
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "31", "dssp": "-", "rsa": 0.531, "bfactor": 18.17, "bfactor_normalized": 0.105, "lip": 1.0, "lip_status": "1", "inter_contacts
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "32", "dssp": "S", "rsa": 0.532, "bfactor": 11.51, "bfactor_normalized": 0.067, "lip": 0.997, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "33", "dssp": "S", "rsa": 0.523, "bfactor": 11.45, "bfactor_normalized": 0.066, "lip": 1.0, "lip_status": "1", "inter_contacts
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "34", "dssp": "-", "rsa": 0.527, "bfactor": 18.39, "bfactor_normalized": 0.106, "lip": 0.96, "lip_status": "1", "inter_contact
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "35", "dssp": "-", "rsa": 0.524, "bfactor": 24.63, "bfactor_normalized": 0.143, "lip": 0.96, "lip_status": "1", "inter_contact
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "36", "dssp": "-", "rsa": 0.516, "bfactor": 31.19, "bfactor_normalized": 0.18, "lip": 0.996, "lip_status": "1", "inter_contact
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "37", "dssp": "-", "rsa": 0.508, "bfactor": 40.38, "bfactor_normalized": 0.234, "lip": 0.892, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "38", "dssp": "H", "rsa": 0.513, "bfactor": 35.36, "bfactor_normalized": 0.205, "lip": 0.996, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "39", "dssp": "H", "rsa": 0.511, "bfactor": 51.6, "bfactor_normalized": 0.299, "lip": 0.888, "lip_status": "1", "inter_contact
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "40", "dssp": "H", "rsa": 0.505, "bfactor": 55.35, "bfactor_normalized": 0.32, "lip": 0.901, "lip_status": "1", "inter_contact
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "41", "dssp": "H", "rsa": 0.492, "bfactor": 41.89, "bfactor_normalized": 0.242, "lip": 0.989, "lip_status": "1", "inter_contac
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "42", "dssp": "H", "rsa": 0.481, "bfactor": 43.2, "bfactor_normalized": 0.25, "lip": 0.993, "lip_status": "1", "inter_contacts
    {"pdb_id": "1jsu", "chain_id": "C", "residue_id": "43", "dssp": "H", "rsa": 0.472, "bfactor": 49.88, "bfactor_normalized": 0.289, "lip": 0.956, "lip_status": "1", "inter_contac
