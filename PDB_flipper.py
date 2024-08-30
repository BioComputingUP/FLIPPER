
import logging

from Bio.PDB import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import FastMMCIFParser

import util_functions as uf


C_RANGE = 3.8
LONG_SHORT_THRESHOLD = 7
MISS_THRESHOLD = 6


# Flipper object to represent a residue
class FlipperResidue:
	# constructor sets important default parameters (has_prev and has_next must be False)
	def __init__(self):
		self.pdb_id = None
		self.model_id = None
		self.pdb_index = None
		self.pdb_insertion_code = None
		self.chain_id = None
		self.pos_in_chain = None
		self.name3 = 'XXX'
		self.name1 = 'X'
		self.has_prev = False
		self.has_next = False
		self.atoms_coord = []
		self.c_alpha_coord = None
		self.uniprot_id = None
		self.uniprot_index = None

	# id of the residue as string with insertion code
	def string_index(self):
		if not self.pdb_insertion_code == ' ':
			return str(self.pdb_index) + self.pdb_insertion_code
		return str(self.pdb_index)

	# !!! RESIDUE FULL IDENTIFIER, used in FLIPPER to identify uniquely a residue. it is used from independent objects as DSSP_MANAGER
	def get_full_identifier(self):
		return "{}_{}_{}_{}{}".format(self.pdb_id, self.model_id, self.chain_id, self.pdb_index, self.pdb_insertion_code)

# Flipper object to represent a chain of residues
# if has a LIST of residues, but also has a DICTIONARY of the residues to access them directly
class FlipperChain:
	# constructor sets the label of the chain and initializes some variables
	def __init__(self, label):
		self.id_label = label
		self.residues = []
		self.string_index_map = {}
		self.rna_dna_chain = False

	# get a residue from his string id (get position in list from map, then retrieve it)
	def get_by_id(self, string_id):
		return self.residues[self.string_index_map[string_id]]

	# get the position in the chain of a residue, given its string id
	def get_pos_by_id(self, string_id):
		return self.string_index_map[string_id]

	# check if a chain contains gaps between two residues.
	def have_gaps(self, res1, res2):
		for i in range(res1.pos_in_chain, res2.pos_in_chain):
			if not self.residues[i].has_next:
				return True
		return False
	

# Flipper object to represent a Structure
# it is a map of labels ==> FlipperChains
class FlipperStructure:
	# constructor
	def __init__(self, pdb_id, mod_id):
		self.pdb_id = pdb_id
		self.model_id = mod_id
		self.chains = {}

	# get all the chains as a list
	def get_chains(self):
		return self.chains.values()


# Network of neighbors in a structure
# it is made as a dictionary, to be used with given methods
class NeighborsNet:
	# constructor
	def __init__(self):
		self.nn = {}

	# add default entry for a residue (empty)
	def add_default(self, res):
		self.nn.setdefault(res.get_full_identifier(), ([],[],[]))

	# add an inter-chain neighbor for res1, that is res2
	def add_inter(self, res1, res2):
		self.nn[res1.get_full_identifier()][0].append(res2) 

	# add a short range intra-chain neighbor for res1, that is res2
	def add_short(self, res1, res2):
		self.nn[res1.get_full_identifier()][1].append(res2) 

	# add a long range intra-chain neighbor for res1, that is res2
	def add_long(self, res1, res2):
		self.nn[res1.get_full_identifier()][2].append(res2) 

	# return all neighbors of a residue
	def get_neighbors(self, residue):
		return self.nn.get(residue.get_full_identifier())

	# return iter chain neighbor for given residue
	def get_inter(self, residue):
		neighbors = self.nn.get(residue.get_full_identifier())
		if neighbors:
			return neighbors[0]
		return None

	# return intra chain neighbor for given residue
	def get_intra(self, residue):
		neighbors = self.nn.get(residue.get_full_identifier())
		if neighbors:
			return neighbors[1]+neighbors[2]
		return None

	# return intra chain short range neighbor for given residue
	def get_short(self, residue):
		neighbors = self.nn.get(residue.get_full_identifier())
		if neighbors:
			return neighbors[1]
		return None

	# return intra chain long range neighbor for given residue
	def get_long(self, residue):
		neighbors = self.nn.get(residue.get_full_identifier())
		if neighbors:
			return neighbors[2]
		return None

	# count inter chain neighbors for given residue
	def get_inter_count(self, residue):
		inter = self.get_inter(residue)
		if inter != None:
			return len(inter)
		return None

	# count intra chain neighbors for given residue
	def get_intra_count(self, residue):
		intra = self.get_intra(residue)
		if intra != None:
			return len(intra)
		return None

	# count short intra chain neighbors for given residue
	def get_short_count(self, residue):
		sh = self.get_short(residue)
		if sh != None:
			return len(sh)
		return None

	# count long intra chain neighbors for given residue
	def get_long_count(self, residue):
		lo = self.get_long(residue)
		if lo != None:
			return len(lo)
		return None

	# count neighbors for given residue
	def get_inter_and_long_count(self, residue):
		neighbors = self.nn.get(residue.get_full_identifier())
		if neighbors != None:
			return len(neighbors[0]), len(neighbors[2])
		return None, None

	# get unique chains from residue neighbors, as pdb and uniprot identifiers
	def get_chains(self, residue):
		chains = []
		uni_chains=[]
		for res in self.get_inter(residue):
			if not (res.chain_id in chains):
				chains.append(res.chain_id)
			if not (res.uniprot_id in uni_chains):
				uni_chains.append(res.uniprot_id)
		return chains, uni_chains


# Structure builder (FROM Bio.PDB objects)
class StructureBuilder:
	# constructor, sets the path of the pdb file and configuration parameters
	def __init__(self, file_path, config):
		self.config = config
		pdb_id = uf.parse_pdb_id_from_file(file_path)
		self.bio_struct = None
		if file_path[-4:] == '.pdb' or file_path[-4:] == '.ent':
			self.bio_struct = PDBParser(QUIET=True).get_structure(pdb_id, file_path)
		elif file_path[-4:] == '.cif':
			self.bio_struct = FastMMCIFParser(auth_chains=False, auth_residues=True, QUIET=True).get_structure(pdb_id, file_path)
			# self.bio_struct = MMCIFParser(QUIET=True).get_structure(pdb_id, file_path)
		if self.bio_struct == None:
			logging.error("{} Flipper unable to parse structure file: {}".format(pdb_id, file_path))
			#sys.exit("ERROR! unable to parse structure: pdb_id: {} file: {}".format(pdb_id, file_path))

	# add a residue to the structure
	def add_residue(self, fl_chain, res):
		# create an empty residue
		fl_res = FlipperResidue()
		### extract res info ###
		fl_res.pdb_id = res_pdb_id(res)
		fl_res.model_id = res_model_id(res)
		fl_res.pdb_index = res_pdb_index(res)
		fl_res.pdb_insertion_code = res_insertion_code(res)
		fl_res.chain_id = fl_chain.id_label
		fl_res.name3 = res.get_resname()
		# if the name in one letter does not exists, let default 'X' 
		n = uf.aa_3to1.get(fl_res.name3)
		if n:
			fl_res.name1 = n
		# if the chain we are inserting in is not a DNA or RNA chain
		if not fl_chain.rna_dna_chain:
			# get the atoms coordinates of the residue, and the alpha carbon also separately
			for atom in res:
				fl_res.atoms_coord.append(atom.get_coord())
				if atom.id == 'CA':
					fl_res.c_alpha_coord = atom.get_coord()
			if fl_res.c_alpha_coord is not None:
				# residues are inserted in order, so give it position equal to the length of the list
				fl_res.pos_in_chain = len(fl_chain.residues)
				# calculate distance between aplah carbons of this new residue and the previous one, to insert gap flags eventually

				if fl_chain.residues and uf.distance_3D(fl_chain.residues[-1].c_alpha_coord, fl_res.c_alpha_coord) < self.config["open_gap_threshold"]:
					fl_chain.residues[-1].has_next = True
					fl_res.has_prev = True
		# if the chain we are inserting in is a DNA or RNA chain
		else:
			# set the uniprot identifier as a string "DNA-RNA"
			fl_res.uniprot_id = "DNA-RNA"
		# add residue inside the chain object
		if fl_res.c_alpha_coord is not None or fl_chain.rna_dna_chain:
			fl_chain.string_index_map[fl_res.string_index()] = len(fl_chain.residues)
			fl_chain.residues.append(fl_res)

	# build the structure, giving it pdb_id as identigier
	def build_structure(self, pdb_id, model_id=0):
		# create Bio.PDB structure
		fl_struct = FlipperStructure(pdb_id, self.bio_struct[model_id].id)

		# for each chain from this structure (first model)
		for chain in self.bio_struct[model_id]:
			# create a FLipper chain with same id
			fl_chain = FlipperChain(chain.id)
			# if this chain is a DNA-RNA chain, set flag
			if is_DNA(chain):
				fl_chain.rna_dna_chain = True
			# for each residue in the chain
			for residue in chain:
				# if it is not hetero, add it (so DNA-RNA residues too)
				if is_good_res(residue):
					self.add_residue(fl_chain, residue)
			# add chain to the structure
			fl_struct.chains[chain.id] = fl_chain
		# return the structure
		return fl_struct

	# create the neighbors network for given struture
	def make_neighbors(self, fl_struct):
		# create an empty NeighborsNet
		nn = NeighborsNet()
		# use NeighborSearch from Bio.PDB to compute distances
		ns = NeighborSearch(list(self.bio_struct.get_atoms()))
		# for each chain in structure that is not a dna-rna one
		for fl_chain in fl_struct.get_chains():
			if not fl_chain.rna_dna_chain:
				# for each residue in this chain
				for fl_res in fl_chain.residues:
					# add a default entry 
					nn.add_default(fl_res)
					# keep track of already inserted nieghbors (the search is mate for each atom in the residue)
					already_have = []
					# for each atom (coordinates) in the residue
					for atom_coord in fl_res.atoms_coord:
						# for each residue in range
						for res in ns.search(atom_coord, self.config["neighbors_range"], level='R'):
							# check if it is good atom and the same model, cause sometimes NS computes all models
							if is_good_res(res) and fl_res.model_id == res_model_id(res) and not res.get_full_id() in already_have:
								# try to get FlipperResidueAssociated
								pos_2 = fl_struct.chains[res_chain_id(res)].string_index_map.get(res_string_index(res))
								# print(fl_res.get_full_identifier(), res.get_full_id(), fl_struct.chains[res_chain_id(res)].string_index_map.get(res_string_index(res)))
								if not pos_2 == None:
									fl_res_2 = fl_struct.chains[res_chain_id(res)].residues[pos_2]
									# if the chain is the same
									if fl_res.chain_id == fl_res_2.chain_id:
										already_have.append(res.get_full_id())
										if fl_res.pos_in_chain == fl_res_2.pos_in_chain:
											continue
										# if distance (as residue number) is less than threshold, then it is a short range neighbor
										if abs(fl_res.pos_in_chain - fl_res_2.pos_in_chain) < self.config["long_short_threshold"] and not fl_chain.have_gaps(fl_res, fl_res_2):
											nn.add_short(fl_res, fl_res_2)
										# else it is a long rage neighbor
										else:
											nn.add_long(fl_res, fl_res_2)
									# if it is not in the same chain it is an inter chain neighbor	
									else:
										nn.add_inter(fl_res, fl_res_2)
										already_have.append(res.get_full_id())

		return nn


# check if a residue is not heteroatom
def is_good_res(res):
	if res.id[0] == ' ':
		return True
	for atom in res:
		if atom.id == 'CA':
			return True
	return False

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

# res string index (id + insertion code) from a Bio.Residue object
def res_string_index(res):
	if not res_insertion_code(res) == ' ':
		return str(res_pdb_index(res)) + res_insertion_code(res)
	return str(res_pdb_index(res))

# get one letter name of an amino acid, Bio.Residue object
def res_name_1_letter(res):
	if res.get_resname() in uf.aa_3to1.keys():
		return uf.aa_3to1[res.get_resname()]
	else:
		return 'X'

# check if given chain is a DNA/RNA chain
def is_DNA(chain):
	# for each residue (cicle it is the only way to access residue info from a chain)
	for res in chain:
		# get the first good residue
		if is_good_res(res):
			return is_dna_res(res)
	# no good residue has been found, discard the chain anyway
	return True

# check if given chain is a DNA/RNA residue
# current strategy is that a DNA/RNA residue is not an hetero atom and does not have an alpha carbon
def is_dna_res(res):
	if is_good_res(res):
		for atom in res:
			if atom.id == 'CA':
				return False
		return True
	return False
	# check if the residue name is one of nucleic acids names (If another different encoding is found add it in this list)
	#return res.get_resname() in [' DA','DU',' DT',' DG',' DC','A','U','T','G','C','  A','  U','  T','  G','  C']
