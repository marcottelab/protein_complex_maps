
#kdrew: this file is for analyzing mass spectrometry data to determine relative stoichiometries between pairs of proteins

import protein_complex_maps.protein_util as pu
import numpy as np
import itertools as it
from scipy.stats import norm
import protein_complex_maps.stoichiometry.stoichiometry as st


#kdrew: function for calculating the ratio of matrix values between two proteins to determine relative stoichiometry
#kdrew: normalized by protein length, takes in uniprot ids
#def calculate_ratio(matrix, name_list, protein_id1, protein_id2, log_transform=True):
def calculate_ratio(msds, protein_id1, protein_id2, log_transform=True):
	#p1_length = pu.get_length_uniprot(protein_id1)
	#p2_length = pu.get_length_uniprot(protein_id2)

	matrix = msds.get_data_matrix()

	array1 = np.array(matrix[msds.get_id_dict()[protein_id1]])[0]
	array2 = np.array(matrix[msds.get_id_dict()[protein_id2]])[0]

	#print array1
	#print array2

	#kdrew: take ratio of matrix values between protein ids
	ratio_list = []
	for a1, a2 in zip(array1, array2):
		#print a1, a2
		if a1 > 0.0 and a2 > 0.0:
			#print a1/a2
			ratio_list.append(a1/a2)

	if log_transform:
		#kdrew: log transform
		ratio_list = np.log10(ratio_list)
		#print ratio_list

	##kdrew: calculate mean and std
	#mean = ratio_list.mean()
	#std = ratio_list.std()
	#print mean, std
    
	#print "ratio_list: ", ratio_list
	return ratio_list


#kdrew: match prot ids to stoichiometry
#kdrew: for all stoichiometries (ex. AB, A2B, AB2, A2B2, A3B, A3B2, A3B3, AB3, A2B3, etc)
#stoichiometry = dict()
#stoichiometry['A'] = 1
#stoichiometry['B'] = 2
#stoichiometry['C'] = 2

#prot_ids = dict()
#prot_ids['A'] = "uniprot_id1"
#prot_ids['B'] = "uniprot_id2"
#prot_ids['C'] = "uniprot_id3"

def relative_stoichiometry_probability( stoichiometry, prior, msds, prot_ids, scale=1.0 ):
	num_data_points = 0
	log_probability = np.log(prior)
	#print "prior log_probability: %s" % (log_probability,)
	#kdrew: for all combinations of pairs (ex. (A,B), (A,C), (B,C))
	for pair in it.combinations(''.join(stoichiometry.keys()), 2):	
 		pair_logratio = np.log10(1.0*stoichiometry[pair[0]]/stoichiometry[pair[1]])
		#print "stoichiometry pair: %s, pair_logratio: %s" % (pair, pair_logratio,)
		pair_norm = norm(loc=pair_logratio, scale=scale)

		ratios = calculate_ratio( msds, prot_ids[pair[0]], prot_ids[pair[1]] )
		print "number of data points between %s and %s : %s" % (prot_ids[pair[0]], prot_ids[pair[1]], len(ratios), )
		num_data_points += len(ratios)
		for r in ratios:
			#print "r: %s : %s : %s" % (r,pair_norm.pdf(r), np.log(pair_norm.pdf(r)),)
			log_probability = log_probability + np.log(pair_norm.pdf(r))
			#print "updated log_probability: %s" % (log_probability,)

	#print stoichiometry, log_probability
	return log_probability, num_data_points



def relative_stoichiometry( msds, ids, stoichiometries ):
	numOfProteins = len(ids)

	stoichiometries_slim = stoichiometries.slim(numOfProteins)

	#print "******"
	#print stoichiometries_slim
	#print "******"

	results = dict()
	num_data_points = 0
	prot_ids = dict()
	#kdrew: assign letter to each protein id
	for i, key in enumerate(stoichiometries_slim[0]):
		prot_ids[key] = ids[i]
	print "prot_ids: %s" % (prot_ids,)

	for stoich in stoichiometries_slim:
		#kdrew: prior is relative to all the other stoichiometries of size numOfProteins (set above)
		prior = stoich.count / stoichiometries_slim.total_count()
		print "stoich: %s prior: %s" % (stoich, prior,)
		log_prob, num_data_points  = relative_stoichiometry_probability( stoich, prior, msds, prot_ids ) 
		results[stoich.__str__()] = log_prob
	
	return results, num_data_points, prot_ids





