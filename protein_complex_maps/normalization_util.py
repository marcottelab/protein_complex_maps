#import math as m
import numpy as np

#kdrew: this function removes rows and columns that are all zero
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
	C = data_matrix + 1.0/M
	return C

#kdrew: adds constant noise over whole matrix, 1.0/M where M = number of rows
def add_noise_over_rows(data_matrix):
	M = data_matrix.shape[0]
	C = data_matrix + 1.0/M
	return C


#kdrew: add in poisson noise to data
#kdrew: lambda parameter is the data_matrix[i,j] value
#kdrew: ability to pass in poisson_module which allows user to set seed, default is unseeded np.random.poisson
def poisson_noise(data_matrix, poisson_module=None):
	poisson_vectorized = None
	if poisson_module == None:
		poisson_vectorized = np.vectorize(np.random.poisson)
	else:
		poisson_vectorized = np.vectorize(poisson_module)

	poisson_mat = poisson_vectorized(data_matrix)

	return poisson_mat

	
