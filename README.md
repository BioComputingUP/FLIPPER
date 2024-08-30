FLIPPER, Fast Linear Interacting Peptides Predictor
==============================================
Damiano Piovesan, Paolo Bonato, Ivan Mičetić and Silvio C.E. Tosatto

Version 2.0

Introduction
------------
This package includes the FLIPPER predictor and the MOBI software to extract LIPs and mobile residues from
PDB structures. FLIPPER is a random forest classifier, MOBI identify mobile residues in NMR ensembles by 
comparing the atomic position in different models of the same PDB. The package is used to generate data in
the **MobiDB database** ([https://mobidb.org](https://mobidb.org))

### References

**MOBI**

Mobi 2.0: An improved method to define intrinsic disorder, mobility and linear binding regions in protein structures \
Piovesan D, Tosatto SCE. \
(2018) Bioinformatics, 34 (1), pp. 122-123 \
[https://pubmed.ncbi.nlm.nih.gov/28968795/](https://pubmed.ncbi.nlm.nih.gov/28968795/)
    
MOBI: A web server to define and visualize structural mobility in NMR protein ensembles \
Martin AJM, Walsh I, Tosatto SCE. \
(2010) Bioinformatics, 26 (22), art. no. btq537, pp. 2916-2917 \
[https://pubmed.ncbi.nlm.nih.gov/20861031/](https://pubmed.ncbi.nlm.nih.gov/20861031/)

**FLIPPER**

FLIPPER: Predicting and Characterizing Linear Interacting Peptides in the Protein Data Bank \
Monzon AM, Bonato P, Necci M, Tosatto SCE, Piovesan D. \
(2021) Jounral of Molecular Biology \
[https://pubmed.ncbi.nlm.nih.gov/33647288/](https://pubmed.ncbi.nlm.nih.gov/33647288/)

Requirements
------------
* Python 3

Usage
-----
The following command prints the help page on screen:

    python3 biodb_disorder.py -h


FLIPPER can be executed by just
providing the PDB input and the output file:

    python3 biodb_disorder.py pdb2zps.ent.gz 2zps.mjson.gz

Configuration
-------------
Check config.ini and config_flipper.json for correct paths to models.
Provide absolute paths. 

Output
------
Output is provided in multiline JSON format, e.g. one JSON document per line. 
Each document corresponds to one residue.

```json
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "25", "dssp": "S", "rsa": 0.501, "bfactor": 84.38, "bfactor_normalized": 0.488, "lip": 0.224, "lip_status": "0", "inter_contacts": 0.909, "intra_long_contacts": 1.455, "helix": 0.0, "beta": 0.455, "coil": 0.545, "delta_rsa": 0.203, "linearity": 0.836, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "26", "dssp": "T", "rsa": 0.502, "bfactor": 86.53, "bfactor_normalized": 0.501, "lip": 0.224, "lip_status": "0", "inter_contacts": 0.818, "intra_long_contacts": 1.545, "helix": 0.0, "beta": 0.455, "coil": 0.545, "delta_rsa": 0.205, "linearity": 0.87, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "27", "dssp": "T", "rsa": 0.499, "bfactor": 77.87, "bfactor_normalized": 0.451, "lip": 0.171, "lip_status": "0", "inter_contacts": 0.818, "intra_long_contacts": 1.636, "helix": 0.0, "beta": 0.455, "coil": 0.545, "delta_rsa": 0.206, "linearity": 0.758, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "28", "dssp": "-", "rsa": 0.486, "bfactor": 51.83, "bfactor_normalized": 0.3, "lip": 0.123, "lip_status": "0", "inter_contacts": 0.818, "intra_long_contacts": 1.727, "helix": 0.0, "beta": 0.455, "coil": 0.545, "delta_rsa": 0.201, "linearity": 0.727, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "29", "dssp": "E", "rsa": 0.486, "bfactor": 24.18, "bfactor_normalized": 0.14, "lip": 0.079, "lip_status": "0", "inter_contacts": 0.636, "intra_long_contacts": 2.0, "helix": 0.0, "beta": 0.545, "coil": 0.455, "delta_rsa": 0.208, "linearity": 0.811, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "30", "dssp": "E", "rsa": 0.479, "bfactor": 21.69, "bfactor_normalized": 0.126, "lip": 0.077, "lip_status": "0", "inter_contacts": 0.545, "intra_long_contacts": 2.182, "helix": 0.0, "beta": 0.636, "coil": 0.364, "delta_rsa": 0.207, "linearity": 0.791, "length_cutoff": 1.0}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "31", "dssp": "E", "rsa": 0.468, "bfactor": 17.82, "bfactor_normalized": 0.103, "lip": 0.054, "lip_status": "0", "inter_contacts": 0.545, "intra_long_contacts": 2.273, "helix": 0.0, "beta": 0.727, "coil": 0.273, "delta_rsa": 0.202, "linearity": 0.699, "length_cutoff": 1.0, "inter_contacts_chains": ["C"]}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "32", "dssp": "E", "rsa": 0.46, "bfactor": 15.93, "bfactor_normalized": 0.092, "lip": 0.054, "lip_status": "0", "inter_contacts": 0.545, "intra_long_contacts": 2.545, "helix": 0.0, "beta": 0.727, "coil": 0.273, "delta_rsa": 0.2, "linearity": 0.718, "length_cutoff": 1.0, "inter_contacts_chains": ["C"]}
{"pdb_id": "1jsu", "chain_id": "A", "residue_id": "33", "dssp": "E", "rsa": 0.457, "bfactor": 23.09, "bfactor_normalized": 0.134, "lip": 0.072, "lip_status": "0", "inter_contacts": 0.545, "intra_long_contacts": 2.545, "helix": 0.0, "beta": 0.818, "coil": 0.182, "delta_rsa": 0.204, "linearity": 0.788, "length_cutoff": 1.0, "inter_contacts_chains": ["C"]}
```

Where fields are:

```text
pdb_id
chain_id
residue_id
dssp
rsa
bfactor
bfactor_normalized
lip
lip_status
inter_contacts
intra_long_contacts
helix
beta
coil
delta_rsa
linearity
length_cutoff
```
