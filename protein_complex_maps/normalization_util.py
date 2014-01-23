
#import math as m
import numpy as np
import math as m
import protein_complex_maps.protein_util as pu

#kdrew: this function removes rows and columns that are all zero
#kdrew: WARNING this removes rows and columns which will invalidate biclusters, only use prior to creating bicluster objects
def remove_zero(data_matrix, zero_rows=True, zero_columns=True):
	clean_data_matrix = data_matrix
	if zero_columns:
		#kdrew: this is a little ugly but np.all produces a truth table where whole column is zero, then negate, then compress (ie. remove those columns)
		clean_data_matrix = clean_data_matrix.compress(~np.array(np.all(clean_data_matrix[:]==0,axis=0))[0],axis=1)
	if zero_rows:
		#kdrew: this is a little ugly but np.all produces a truth table where whole row is zero, then negate, then compress  (i.e remove those rows)
		clean_data_matrix = clean_data_matrix.compress(~np.array(np.all(clean_data_matrix[:]==0,axis=1)).reshape(-1),axis=0)

	return clean_data_matrix

#kdrew: converts cells to frequencies where whole row sums to one
def normalize_over_columns(data_matrix):
	return np.nan_to_num( data_matrix / np.sum(data_matrix,1) )

#kdrew: converts cells to frequencies where whole column sums to one
def normalize_over_rows(data_matrix):
	return np.nan_to_num( data_matrix / np.sum(data_matrix,0) )


#kdrew: adds constant noise over whole matrix, 1.0/M where M = number of columns
def add_noise_over_columns(data_matrix):
	M = data_matrix.shape[1]
	C = add_noise(data_matrix, 1.0/M )
	return C

#kdrew: adds constant noise over whole matrix, 1.0/M where M = number of rows
def add_noise_over_rows(data_matrix):
	M = data_matrix.shape[0]
	C = add_noise(data_matrix, 1.0/M )
	return C

def add_noise(data_matrix, noise_constant):
	C = data_matrix + noise_constant
	return C


#kdrew: add in noise to data
#kdrew: default poisson and lambda parameter is the data_matrix[i,j] value
#kdrew: ability to pass in sample_module which allows user to set seed and sampling model, default is unseeded np.random.poisson
def sample_noise(data_matrix, sample_module=None):
	sample_vectorized = None
	if sample_module == None:
		sample_vectorized = np.vectorize(np.random.poisson)
	else:
		sample_vectorized = np.vectorize(sample_module)

	sample_mat = sample_vectorized(data_matrix)

	return sample_mat


def log_transform(data_matrix, base=m.e):
	log_vec = np.vectorize(m.log)
	return log_vec(data_matrix, base)
	
def exp_transform(data_matrix):
	exp_vec = np.vectorize(m.exp)
	return exp_vec(data_matrix)
	
#kdrew: function scales all values to be above 1.0
#kdrew: find min value and multiplies all values 1/min
def min_to_one_scale(data_matrix):
	min = np.min(data_matrix)
	return data_matrix/(1.0*min)

#kdrew: if values are below threshold, entry is marked 0
def threshold(data_matrix, threshold, absent=0):
	below_threshold = data_matrix < threshold

	data_matrix[below_threshold] = absent

	print "number below_threshold %s" % (len(below_threshold),)

	return data_matrix


#kdrew: if values are above threshold, entry is marked present (1) otherwise absent (0)
def binary(data_matrix, threshold, present=1, absent=0):
	above_threshold = data_matrix >= threshold
	below_threshold = data_matrix < threshold

	data_matrix[above_threshold] = present
	data_matrix[below_threshold] = absent

	return data_matrix

#kdrew: normalizes by protein length, if no length available make NAN
def normalize_length(data_matrix, id_dict, initialize=np.nan):
	length_array = np.repeat(initialize, data_matrix.shape[0])

	length_dict = pu.get_length_uniprot(id_dict.keys())
	for id in id_dict.keys():
		try:
			length_array[id_dict[id]] = length_dict[id]
		except KeyError:
			continue

	print data_matrix
	#kdrew: transpose to do division
	data_matrix = data_matrix.T/length_array

	print data_matrix.T
	#kdrew: "retranspose" to get data matrix back into form
	return data_matrix.T

	
