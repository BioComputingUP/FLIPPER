#!/usr/bin/env python

import logging
import os
import pickle

import numpy as np

import PDB_flipper
import util_functions as uf


# Class that manages the lips_dataset.txt file, from parsing to features extraction
class LipsFinder:
	# constructor, set configurations, build a dssp manager object
	def __init__(self, configurations, dssp_dict, dssp_dict_chains):
		self.config = configurations
		self.dssp = dssp_dict  # {<residue_full_id>: (<ss>, <asa>)}
		self.chain_dssp = dssp_dict_chains  # same as dssp but caclulated on splitted chains
		self.lip_classifier = None
		#print(self.dssp)

	# extract features from a structure, using neighbors and dssp 
	# for training pourpose the output can be flattened twice
	def extract_features(self, structure, neighbors, pdb_id, file_path, with_lab=False, dataset=None):
		logging.debug("{} - starting features extraction".format(pdb_id))
		# check that with lab is used properly
		if with_lab and dataset == None:
			logging.debug("{} - requested truth labels extraction but no dataset is given".format(pdb_id))

		# select the maximum window size
		max_win = max(self.config["asa_ss_window"],self.config["neighbors_window"],self.config["linearity_window"])
		
		# init dict of chains extracted
		extracted_pdb_X = {}
		extracted_pdb_y = {}
		extracted_pdb_r = {}

		# for each chain in the structure
		for chain in structure.get_chains():
			# if with lab flag is enabled, skip unlabelled chains
			if (not with_lab or chain.id_label in dataset.get_labelled_chains(pdb_id)) and not chain.rna_dna_chain:

				# init list of residues extracted
				extracted_pdb_X[chain.id_label] = np.zeros(shape=(len(chain.residues), 11))
				extracted_pdb_y[chain.id_label] = np.empty(len(chain.residues), dtype=int)

				# for each of the residues in the chain
				for res in chain.residues:

					# init value of extraction for current residue
					extracted_residue_X = np.zeros(11)
					extracted_residue_y = 0

					###### EXTEND THE WINDOW ######
					# big as 2*max_win + 1 or stops when encountering gaps
					res_index = res.pos_in_chain
					start = res_index
					stop = res_index + 1
					while start > res_index - max_win and start > 0 and chain.residues[res_index].has_prev:
						start -= 1
					while stop < res_index + max_win + 1 and stop < len(chain.residues) and chain.residues[res_index].has_next:
						stop += 1

					###### NEIGHBORS & SS WINDOW ######
					# limit max window to feature specific range
					neigh_start = max((res_index - self.config["neighbors_window"]), start)
					neigh_stop = min((res_index + self.config["neighbors_window"] + 1), stop)

					for r in chain.residues[neigh_start : neigh_stop]:
						inter, intra_long = neighbors.get_inter_and_long_count(r)
						extracted_residue_X[0] += inter
						extracted_residue_X[1] += intra_long

					# add inter chain / long intra chain
					extracted_residue_X[0] = extracted_residue_X[0]/(neigh_stop - neigh_start)
					extracted_residue_X[1] = extracted_residue_X[1]/(neigh_stop - neigh_start)

					count = 0
					for r in chain.residues[neigh_start: neigh_stop]:
						# _tmp = r.get_full_identifier()
						entry = self.dssp.get(r.get_full_identifier())
						if entry:
							# get data from dssp
							ss = entry[0]
							#add helices
							if ss in ['H', 'G', 'I']:
								extracted_residue_X[2] += 1
							# add beta sheets
							elif ss in ['E', 'B']:
								extracted_residue_X[3] += 1
							# add non ss
							else:
								extracted_residue_X[4] += 1
							count += 1
					if count > 0:
						extracted_residue_X[2] = extracted_residue_X[2]/(count)
						extracted_residue_X[3] = extracted_residue_X[3]/(count)
						extracted_residue_X[4] = extracted_residue_X[4]/(count)

					###### ASA WINDOW ######
					# limit max window to feature specific range
					asa_start = max((res_index - self.config["asa_ss_window"]), start)
					asa_stop = min((res_index + self.config["asa_ss_window"] + 1), stop)

					count = 0
					for r in chain.residues[asa_start: asa_stop]:

						pdb_entry = self.dssp.get(r.get_full_identifier())
						chain_entry = self.chain_dssp.get(r.get_full_identifier())
						# print(r.get_full_identifier(), pdb_entry, chain_entry)
						if pdb_entry and chain_entry:
							# get data from dssp
							asa = pdb_entry[1]
							chain_asa = chain_entry[1]
							# asa can be a string 'NA' with DNA 
							if not chain_asa == 'NA' and not asa == 'NA':
								extracted_residue_X[5] += chain_asa
								extracted_residue_X[6] += (chain_asa - asa)
							count += 1

					# check number of computed residues in window
					if count > 0:
						extracted_residue_X[5] = extracted_residue_X[5]/(count)
						extracted_residue_X[6] = extracted_residue_X[6]/(count)

					###### MEDIUM WINDOW AND OTHER ######
					# limit max window to feature specific range
					lin_start = max((res_index - self.config["linearity_window"]), start)
					lin_stop = min((res_index + self.config["linearity_window"] + 1), stop)

					extracted_residue_X[7] = neighbors.get_inter_count(res)
					extracted_residue_X[8] = neighbors.get_long_count(res)
					extracted_residue_X[9] = uf.distance_3D(chain.residues[lin_start].c_alpha_coord, chain.residues[lin_stop-1].c_alpha_coord)/(lin_stop-lin_start)
					extracted_residue_X[10] = min(self.config["chain_length_limit"], len(chain.residues))/self.config["chain_length_limit"]
					# if with lab falg is enabled check if residue is lip or not
					if with_lab and dataset.check_res(res):
						extracted_residue_y = 1
					# append extracted residue to residues list
					extracted_pdb_X[chain.id_label][res_index] = extracted_residue_X
					extracted_pdb_y[chain.id_label][res_index] = extracted_residue_y

		# return extraction
		"""
		Features (extracted_residue_X)
		0 inter contacts (window neighbors_window)
		1 intra long range contacts (window neighbors_window)
		2 helix (window neighbors_window)
		3 beta (window neighbors_window)
		4 non-ss (window neighbors_window)
		5 rsa (window asa_ss_window)
		6 delta-rsa (window asa_ss_window)
		7 inter contacts
		8 long range contacts
		9 distance_3D (window linearity_window)
		10 length cutoff (cap chain_length_limit)
		"""
		return extracted_pdb_X, extracted_pdb_y

	# output the trained model to file
	def model_out(self, path = None):
		if path == None:
			path = os.path.join(os.path.dirname(__file__), self.config["trained_model_path"])
		logging.debug("dumping trained model to file: {}".format(path))
		pickle.dump(self.lip_classifier, open(path, 'wb'))
		return 

	# parse trained model from file
	def model_in(self, path = None):
		if path == None:
			path = os.path.join(os.path.dirname(__file__), self.config["trained_model_path"])
		logging.debug("loading trained model from file: {}".format(path))
		self.lip_classifier = pickle.load(open(path, 'rb'))
		return

	# predict a pdb file
	def predict(self, pdb_id, file_path, proba = False, blur = 0, threshold = 0.5,  gap = 0, model_id = 0):
		logging.debug("{} - prediction started from file: {}".format(pdb_id, file_path))
		# init results
		raw_predictions = {}
		predictions = {}

		# build a StructureBuilder obj with configurations, then obtain structure and neighbors
		builder = PDB_flipper.StructureBuilder(file_path, self.config)

		structure = builder.build_structure(pdb_id, model_id)
		neigh = builder.make_neighbors(structure)

		# if there is only one chain exit prediction
		if len(structure.get_chains()) < 2:
			logging.debug("{} - only one chain found, skipping prediction.".format(pdb_id))
			#sys.exit("INFO: only one chain in pdb {}. skipping prediction.".format(pdb_id))
			return predictions, raw_predictions, structure, neigh, None

		# check if there are only dna/rna chains. if true, exit prediction
		only_dna_rna = True
		for c in structure.get_chains():
			if not c.rna_dna_chain:
				only_dna_rna = False
		if only_dna_rna:
			logging.debug("{} - only DNA-RNA chains found, skipping prediction.".format(pdb_id))
			#sys.exit("INFO: only one chain in pdb {}. skipping prediction.".format(pdb_id))
			return predictions, raw_predictions, structure, neigh, None

		# extract features from the pdb
		X, y = self.extract_features(structure, neigh, pdb_id, file_path, with_lab = False)
		# logging.debug("{} - features extracted".format(pdb_id))
		# logging.debug("X: {}, y: {}".format(X, y))
		# for each chain, perform post processing operations
		# logging.debug("classifier {}".format(self.lip_classifier))
		for chain_id in X.keys():
			if not proba:
				# logging.debug("{} - predicting: {}_{}".format(pdb_id, pdb_id, chain_id))
				raw_predictions[chain_id] = self.lip_classifier.predict(X[chain_id])
				# logging.debug("{} - post_processing: {}_{}".format(pdb_id, pdb_id, chain_id))
				predictions[chain_id] = PDB_flipper.gap_fill(raw_predictions[chain_id], g=gap)
			else:
				# logging.debug("{} - predicting: {}_{}".format(pdb_id, pdb_id,chain_id))
				raw_predictions[chain_id] = uf.blur_predictions(self.lip_classifier.predict_proba(X[chain_id])[:,1], w=blur)
				# logging.debug("{} - post_processing: {}_{}".format(pdb_id, pdb_id,chain_id))
				predictions[chain_id] = uf.gap_fill(np.array([1 if res_score >= threshold else 0 for res_score in raw_predictions[chain_id]]), g=gap)
		# return processed predictions, original predictions, structure and neighbors network
		return predictions, raw_predictions, structure, neigh, X


