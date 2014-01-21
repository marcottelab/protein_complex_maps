
#kdrew: this file is for analyzing mass spectrometry data to determine relative stoichiometries between pairs of proteins

import protein_complex_maps.protein_util as pu
import numpy as np
import itertools as it
from scipy.stats import norm
import protein_complex_maps.stoichiometry.stoichiometry as st


#stoichiometries = [[1,1],[2,1],[2,2],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[4,4],[5,1],[5,5],[6,1],[6,2],[6,3],[6,6],[7,7],[8,1],[8,2],[8,4],[8,6],[8,8],[9,1],[9,3],[9,9],[10,5],[10,10],[12,1],[12,2],[12,6],[12,12],[14,7],[14,14],[15,12],[16,8],[18,12],[20,20],[24,12],[24,24],[60,48],[60,60],[64,32],[108,108],[120,60],[132,33],[144,36],[168,42],[180,60],[180,120],[180,180],[240,240],[360,180],[420,420],[780,60],[780,120],[1,1,1],[2,1,1],[2,2,1],[2,2,2],[3,1,1],[3,3,1],[3,3,2],[3,3,3],[4,1,1],[4,2,2],[4,3,1],[4,4,1],[4,4,4],[5,1,1],[5,2,1],[5,3,2],[5,5,5],[6,3,3],[6,4,2],[6,4,4],[6,6,3],[8,6,6],[10,10,1],[10,10,10],[12,2,2],[14,14,14],[18,18,12],[24,24,24],[36,12,6],[60,60,60],[120,120,60],[180,120,120],[240,180,180],[240,240,240],[480,240,240],[720,120,60],[780,120,60],[900,60,60],[2,1,1,1],[2,2,1,1],[2,2,2,1],[2,2,2,2],[3,1,1,1],[3,2,1,1],[3,3,1,1],[3,3,2,2],[3,3,3,3],[4,2,2,2],[4,4,1,1],[4,4,4,2],[4,4,4,4],[6,3,3,3],[8,8,8,8],[12,6,6,6],[14,4,4,4],[14,5,5,5],[14,6,5,5],[14,6,6,6],[18,18,6,3],[60,60,60,60],[72,72,36,36],[120,60,60,60],[176,176,176,176],[180,180,120,120],[180,180,180,180],[240,60,60,60],[600,240,60,60],[720,60,60,60],[1200,120,120,60],[2,1,1,1,1],[2,2,2,2,1],[2,2,2,2,2],[3,3,1,1,1],[7,2,2,2,1],[10,3,3,1,1],[10,8,4,4,2],[12,6,6,6,6],[16,8,4,4,4],[60,60,60,60,60],[120,60,60,60,60],[180,60,60,60,60],[240,240,240,240,240],[720,240,120,60,60],[2,1,1,1,1,1],[2,2,1,1,1,1],[3,1,1,1,1,1],[3,3,1,1,1,1],[8,3,3,1,1,1],[10,3,3,1,1,1],[60,60,60,60,60,60],[600,600,120,60,60,60],[600,600,180,120,60,60],[3,3,2,2,1,1,1],[4,2,2,2,2,2,2],[2,2,2,2,2,2,2,1],[2,2,2,2,2,2,2,2],[4,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2],[3,3,1,1,1,1,1,1,1],[4,2,2,2,2,2,2,2,2],[4,2,2,2,2,2,2,2,2,2],[10,6,4,2,2,2,2,2,2,2],[20,6,6,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,1,1],[6,3,3,3,3,3,3,3,3,3,3],[1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,1],[1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,4,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,4,4,2,2,2,2,2,2,2,2,2,2,2],[14,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

#kdrew: function for calculating the ratio of matrix values between two proteins to determine relative stoichiometry
#kdrew: normalized by protein length, takes in uniprot ids
#def calculate_ratio(matrix, name_list, protein_id1, protein_id2, log_transform=True):
def calculate_ratio(msds, protein_id1, protein_id2, log_transform=True):
	#p1_length = pu.get_length_uniprot(protein_id1)
	#p2_length = pu.get_length_uniprot(protein_id2)

	matrix = msds.get_data_matrix()

	array1 = np.array(matrix[msds.get_id_dict()[protein_id1]])[0]
	array2 = np.array(matrix[msds.get_id_dict()[protein_id2]])[0]

	print array1
	print array2

	#kdrew: take ratio of matrix values between protein ids
	ratio_list = []
	for a1, a2 in zip(array1, array2):
		print a1, a2
		if a1 > 0.0 and a2 > 0.0:
			print a1/a2
			ratio_list.append(a1/a2)

	if log_transform:
		#kdrew: log transform
		ratio_list = np.log10(ratio_list)
		print ratio_list

	##kdrew: calculate mean and std
	#mean = ratio_list.mean()
	#std = ratio_list.std()
	#print mean, std
    
	print "ratio_list: ", ratio_list
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
	log_probability = np.log(prior)
	print "log_probability: %s" % (log_probability,)
	#kdrew: for all combinations of pairs (ex. (A,B), (A,C), (B,C))
	for pair in it.combinations(''.join(stoichiometry.keys()), 2):	
		print pair
 		pair_logratio = np.log10(1.0*stoichiometry[pair[0]]/stoichiometry[pair[1]])
		print "pair_logratio: %s" % (pair_logratio, )
		pair_norm = norm(loc=pair_logratio, scale=scale)

		ratios = calculate_ratio( msds, prot_ids[pair[0]], prot_ids[pair[1]] )
		for r in ratios:
			print "r: %s : %s : %s" % (r,pair_norm.pdf(r), np.log(pair_norm.pdf(r)),)
			log_probability = log_probability + np.log(pair_norm.pdf(r))
			print "log_probability: %s" % (log_probability,)

	return log_probability



def relative_stoichiometry( msds, ids, stoichiometries ):
	numOfProteins = len(ids)

	l=[ x for x in stoichiometries if len(x) == numOfProteins]
	print l
	stoichiometries_slim = st.Stoichiometries(l=l)

	print "******"
	print stoichiometries_slim

	results = dict()
	for stoich in stoichiometries_slim:
		prot_ids = dict()
		#kdrew: assign letter to each protein id
		for i, key in enumerate(stoich):
			prot_ids[key] = ids[i]

		#kdrew: prior is relative to all the other stoichiometries of size numOfProteins (set above)
		prior = stoich.count / stoichiometries_slim.total_count()
		log_prob = relative_stoichiometry_probability( stoich, prior, msds, prot_ids ) 
		results[stoich.__str__()] = log_prob
	
	return results





