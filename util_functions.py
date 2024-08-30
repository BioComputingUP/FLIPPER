#!/usr/bin/env python

import re
import gzip
import numpy as np
import xml.etree.ElementTree as ET
import logging


aa_3to1 = {
	'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
	#'GLX': 'Z', 'ASX': 'B', 'TER': '*', 'XAA': 'X'
}

aa_1to3 = {
	'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
	'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
	'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
	'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'
	#'Z': 'GLX', 'B': 'ASX', '*': 'TER', 'X': 'XAA'
}


#######################################################
################  Residues Help Func  #################
#######################################################

# check if a residue is not heteroatom
def is_good_res(res):
	return res.id[0] == ' '

# extract the index of the residue (relative to pdb file)
def res_pdb_index(res):
	return res.id[1]

# extract the insertion of the residue
def res_insertion_code(res):
	return res.id[2]

# extract the id of the pdb from a residue
def res_pdb_id(res):
	return res.get_full_id()[0]

# extract the id of the model from a residue
def res_model_id(res):
	return res.get_full_id()[1]

# extract the id of the chain from a residue
def res_chain_id(res):
	return res.get_full_id()[2]

def res_string_index(res):
	if not res_insertion_code(res) == ' ':
		return str(res_pdb_index(res)) + res_insertion_code(res)
	return str(res_pdb_index(res))

def res_name_1_letter(res):
	if res.get_resname() in aa_3to1.keys():
		return aa_3to1[res.get_resname()]
	else:
		return 'X'

# check if given chain is a DNA chain
# (in some pdb files they are very long to compute and DSSP would print warnings)
def is_DNA(chain):
	# for each residue (cicle it is the only way to access residue info from a chain)
	for res in chain:
		# get the first good residue
		if is_good_res(res):
			return is_dna_res(res)
	# no good residue has been found, discard the chain anyway
	return True

def is_dna_res(res):
	# check if the residue name is one of nucleic acids names (If another different encoding is found add it in this list)
	return res.get_resname() in [' DA','DU',' DT',' DG',' DC','A','U','T','G','C','  A','  U','  T','  G','  C']


def parse_missing_from_header(pdb_path):
	file = open(pdb_path, 'r')
	lines = file.read().splitlines()
	missing = {}
	found = False
	for line in lines:
		if re.match("^REMARK\s*465\s*[A-Z]{3}\s*\w\s*\d+\s*", line):
			found = True
			tokens = line.split()
			missing.setdefault(tokens[3], [])
			if re.match("^\d*$",tokens[4]):
				missing[tokens[3]].append(int(tokens[4]))
			else:
				missing[tokens[3]].append(int(tokens[4][0:-1]))
		elif found == True or re.match("^ATOM.*", line):
			break
	file.close()
	return missing


# get the spatial center of the residue (mean between all atoms)
def get_center(residue):
		cm = []
		# collect coordinates of every atom in residue
		for atom in residue:
			cm.append(atom.get_coord())
		# sum coord and divide for number of atoms
		return np.sum(np.array(cm), axis=0) / len(cm)

def get_chain_center(residues):
	cm = []
	for residue in residues:
		cm.append(get_center(residue))
	return np.sum(np.array(cm), axis=0)/len(cm)

def get_chain_radius_error(residues):
	center = get_chain_center(residues)
	radius_list = []
	for residue in residues:
		radius_list.append(np.linalg.norm(center - get_center(residue)))
	mean_radius = np.sum(np.array(radius_list), axis=0)/len(radius_list)
	errors = []
	for residue in residues:
		errors.append((mean_radius-np.linalg.norm(center - get_center(residue)))/mean_radius)
	return errors

def distance_3D(point_1, point_2):
	return np.linalg.norm(point_1-point_2)

# compute distance between two residues
def compute_distance(res1, res2):
	# initialize coordinates variables
	a = 0
	b = 0

	# search for the Carbon Alpha atom and get coordinates
	for atom in res1:
		if atom.id == "CA":
			a = atom.get_coord()

	for atom in res2:
		if atom.id == "CA":
			b = atom.get_coord()

	# return norm as the distance
	return distance_3D(a, b)


# compute angle between three residues
def compute_angle(res1, res2, res3):
	# initialize coordinates variables
	a = 0
	b = 0
	c = 0

	# search for the Carbon Alpha atom and get coordinates
	for atom in res1:
		if atom.id == "CA":
			a = atom.get_coord()
	for atom in res2:
		if atom.id == "CA":
			b = atom.get_coord()
	for atom in res3:
		if atom.id == "CA":
			c = atom.get_coord()

	# compute the edges ba and bc
	ba = a - b
	bc = c - b
	# with a very small window can happens that one of the edges is 0
	if (np.linalg.norm(ba) * np.linalg.norm(bc) == 0):
		# print a warning and return a 90 degrees angle
		print("Wrong arguments calculating angle between: {}, {}, {}".format(res1.get_full_id(),res2.get_full_id(),res3.get_full_id()))
		return 90
	# if not 0 compute and return angle
	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle = np.arccos(cosine_angle)
	return np.degrees(angle)

# fill gaps in predictions g is the maximum number of residues to consider a region as a gap
def gap_fill(chain_pred, g=0):
	# if g = 0 retur the predictions as they are
	if(g == 0):
		return chain_pred
	# start with the positives
	for index in range(len(chain_pred)):
		# if the residue in index position is positive
		if chain_pred[index] == 1:
			# init count with one residue
			count = 1
			#start from residue to the left
			i = index - 1
			# while the index is inside the chain and the prediction is 1 and the count is less or equal to g
			while i >= 0 and chain_pred[i] == 1 and count <= g:
				# increment count and move left
				count += 1
				i -= 1
			# then go right
			i = index + 1
			# while the index is inside the chain and the prediction is 1 and the count is less or equal to g (can be greater already)
			while i < len(chain_pred) and chain_pred[i] == 1 and count <= g:
				# increment count and move right
				count += 1
				i += 1	
			# if the count is less or equal to g, it is a gap	
			if count <= g:
				chain_pred[index] = 0
	# for negatives is a little bit different (if at beginning or end it is not a gap)
	for index in range(len(chain_pred)):
		# if the residue in index position is negative
		if chain_pred[index] == 0:
			# init count with one residue
			count = 1
			#start from residue to the left
			i = index - 1
			# while the index is inside the chain and the prediction is 0 and the count is less or equal to g
			while i >= 0 and chain_pred[i] == 0 and count <= g:
				# increment count and move left
				count += 1
				i -= 1
			# if the beginning has been reach i should be equal to -1. then leave the 0 as it is
			if i < 0:
				continue
			# then go right
			i = index + 1
			# while the index is inside the chain and the prediction is 0 and the count is less or equal to g
			while i < len(chain_pred) and chain_pred[i] == 0 and count <= g:
				# increment count and move right
				count += 1
				i += 1
			# if the end has been reach i should be equal to len(chain_pred). then leave the 0 as it is
			if i >= len(chain_pred):
				continue
			# if beginning and end not reached and g is less or equal to gap number, set prediction to 1		
			if count <= g:
				chain_pred[index] = 1
	# return the gap_filled predictions 
	return chain_pred


def blur_predictions(chain_pred, w=0):
	if(w == 0):
		return chain_pred
	new_pred = []
	for index in range(len(chain_pred)):
		start = max(0, index - w)
		stop  = min(len(chain_pred), index + w + 1)
		new_pred.append(np.sum(chain_pred[start:stop])/(stop-start))
	return np.array(new_pred)


def find_regions(labels):
	regions = []
	start = None
	for l_idx in range(len(labels)):
		if labels[l_idx] == '1' and start == None:
			start = l_idx
		if labels[l_idx] == '0' and start != None:
			regions.append((start, l_idx-1))
			start = None
	if start != None:
		regions.append((start, len(labels)-1))
	return regions


def parse_mapping(file_path):
	mapping = {}
	tree = ET.parse(file_path)
	root = tree.getroot()
	
	# all item attributes
	for elem in root:
		if 'type' in elem.attrib and elem.attrib['type'] == "protein":
			for segment in elem:
				for residue in segment[0]:
					pdb_idx = None
					chain_id = None
					uniprot_idx = None
					uniprot_id = None
					for entry in residue:
						if entry.attrib["dbSource"] == "PDB" and not entry.attrib["dbResNum"] == "null":
							pdb_idx = entry.attrib["dbResNum"]
							chain_id = entry.attrib['dbChainId']
						elif entry.attrib["dbSource"] == "UniProt":
							uniprot_id = entry.attrib["dbAccessionId"]
							uniprot_idx = entry.attrib["dbResNum"] 
					if chain_id and pdb_idx and uniprot_idx and uniprot_id:
						mapping.setdefault(chain_id, {})
						mapping[chain_id].setdefault(uniprot_id, {})
						mapping[chain_id][uniprot_id][pdb_idx] = uniprot_idx
	return mapping


def map_structure(structure, mapping):
	for chain_id in mapping.keys():
		chain = structure.chains.get(chain_id)
		if chain:
			for uniprot_id in mapping[chain_id].keys():
				for pdb_idx in mapping[chain_id][uniprot_id]:
					res = chain.get_by_id(pdb_idx)
					if res:
						res.uniprot_id = uniprot_id
						res.uniprot_index = mapping[chain_id][uniprot_id][pdb_idx]
	return mapping


def is_gz(file_path):
	return file_path[-3:] == '.gz'

def get_suffix(file_path):
	if is_gz(file_path):
		return file_path[-7:-3]
	else:
		return file_path[-4:]
	

def extract_gz(path_from, path_to):
	inputfile = gzip.GzipFile(path_from, 'rb')
	s = inputfile.read()
	inputfile.close()
	outputfile = open(path_to, 'wb')
	outputfile.write(s)
	outputfile.close()
	return

def make_gz(path_from, path_to):
	inputfile = open(path_from, 'rb')
	s = inputfile.read()
	inputfile.close()
	outputfile = gzip.GzipFile(path_to, 'wb')
	outputfile.write(s)
	outputfile.close()
	return


def parse_pdb_id_from_file(file_path):
	if '.ent' in file_path or '.pdb' in file_path:
		file = open(file_path)
		lines = file.read().splitlines()
		for line in lines:
			tokens = line.split()
			if tokens[0] == "HEADER":
				return tokens[-1].lower()
	elif '.cif' in file_path:
		file = open(file_path)
		lines = file.read().splitlines()
		for line in lines:
			tokens = line.split()
			if tokens[0] == "_entry.id":
				return tokens[-1].lower()
	else:
		logging.warning("Can't extract PDB id from file " + file_path)
		return 'noid'

def get_ranges(labels):
	ranges = []
	start = None
	for idx in range(len(labels)):
		if not start == None:
			if labels[idx] == 0:
				ranges.append((start, idx-1))
				start = None
		else:
			if labels[idx] == 1:
				start = idx
	if not start == None:
		ranges.append((start, len(labels)-1))
	return ranges


