
import numpy as np
import protein_complex_maps.normalization_util as nu

def multiple_dot(matrix):

	if matrix.shape[1] == 0 or matrix.shape[0] == 0:
		return 0.0

	result_vector = np.ones(matrix.shape[1])
	for row in matrix:
		row = np.asarray(row).reshape(-1)
		#print "row: %s" % (row,)
		#print "before result_vector: %s" % (result_vector,)
		result_vector = result_vector * row
		#print "result_vector: %s" % (result_vector,)

	total_sum = np.sum(result_vector)
	#print "total_sum: %s" % (total_sum,)

	#return total_sum/(matrix.shape[1] * matrix.shape[0])
	return total_sum

def multiple_dot_per_unit(matrix):
	total_score = multiple_dot(matrix)
	return total_score/(matrix.shape[1] * matrix.shape[0])


def multiple_dot_neg(matrix):
	return -1.0*multiple_dot(matrix)

def sum_cells(matrix):
	#kdrew: if there are no rows or no columns, return large negative number (i.e. bad score)
	#kdrew: should I raise an exception here?
	if matrix.shape[1] == 0 or matrix.shape[0] == 0:
		return -999999999.0

	log_dmat = nu.log_transform(matrix)
	#print log_dmat

	log_sum = sum(log_dmat)
	#print "log_sum: %s" % (log_sum,)
	exp_sum = nu.exp_transform(log_sum)
	#print "exp_sum: %s" % (exp_sum,)

	#total_sum = sum(sum(np.array(log_dmat.reshape(-1))))
	total_sum = sum(np.array(exp_sum.reshape(-1)))

	#print "total_sum: %s" % (total_sum,)
	return total_sum/(matrix.shape[1] * matrix.shape[0])


#kdrew: deprecated, use random_sampling_util module directly
def zscore(matrix, full_data_matrix, score_function, randsamp):
	#kdrew: get the mean and std of random biclusters the same size as input matrix
	mean, std = randsamp.score_mean_std_all( full_data_matrix, matrix.shape[0], matrix.shape[1] )
	score = score_function(matrix)
	zscore = abs(mean - score)/std

	return zscore

