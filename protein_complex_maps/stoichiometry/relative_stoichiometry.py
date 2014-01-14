
#kdrew: this file is for analyzing mass spectrometry data to determine relative stoichiometries between pairs of proteins

import protein_complex_maps.protein_util as pu
import numpy as np

#kdrew: function for calculating the ratio of matrix values between two proteins to determine relative stoichiometry
#kdrew: normalized by protein length, takes in uniprot ids
def calculate_ratio(matrix, name_list, protein_id1, protein_id2, log_transform=True):
	p1_length = pu.get_length_uniprot(protein_id1)
	p2_length = pu.get_length_uniprot(protein_id2)

	array1 = np.array(matrix[name_list.index(protein_id1)])[0]
	array2 = np.array(matrix[name_list.index(protein_id2)])[0]

	print array1
	print array2

	#kdrew: normalize matrix with respect to length
	array1 = array1/p1_length
	array2 = array2/p2_length

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






