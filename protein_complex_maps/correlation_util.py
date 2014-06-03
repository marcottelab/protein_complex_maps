from scipy.stats.stats import pearsonr
import math as m
import numpy as np
import protein_complex_maps.normalization_util as nu
import argparse
import pickle



CLOSE_TO_ONE = 0.99999999
EPSILON = 0.000001

def main():

	parser = argparse.ArgumentParser(description="Hierarchical clusters fractionation data")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--sample_method", action="store", dest="sample_method", required=False, default=None,
						help="Sampling method to add noise to correlation calculations (poisson or normal)")
	parser.add_argument("--sampling_iterations", action="store", type=int, dest="average_cnt", required=False, default=3,
						help="Number of samples to compute average correlation")
	parser.add_argument("--one_vs_all", action="store_true", dest="one_vs_all", required=False, default=False,
						help="Compares a single protein to all other proteins in msds")

	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )

	if args.sample_method == "poisson":
		sample_module = np.random.poisson
	elif args.sample_method == "normal":
		sample_module = np.random.normal 
	else:
		sample_module = None
	
	if args.one_vs_all:
		#kdrew: take first protein
		index = msds.get_id_dict()[args.proteins[0]]
		corr_list = correlation_array(msds.get_data_matrix(), index)
		print corr_list



	for prot in args.proteins:
		index = msds.get_id_dict()[prot]

		scores, tvals = sample_correlation_distribution(matrix=msds.get_data_matrix(), index=index, iterations=args.average_cnt, sample_module=sample_module)

		print scores

def correlation_array (matrix, index):
	corr_list = []

	data_array1 = np.array(matrix[index].reshape(-1))[0]
	for i in range(matrix.shape[0]):
		data_array2 = np.array(matrix[i].reshape(-1))[0]
		corr_list.append(pearsonr(data_array1, data_array2))

	return corr_list
		

def correlation_distribution( matrix, index=None ):

	corrcoefMat = np.nan_to_num(np.corrcoef(matrix))
	#print corrcoefMat
	if index == None:
		#kdrew: gets indices of upper triangle (w/o diagonal) and returns values in a list
		#score_distribution = corrcoefMat[np.triu_indices(len(corrcoefMat),1)].tolist()
		score_distribution = corrcoefMat[np.triu_indices(len(corrcoefMat),1)]

	else:
		#kdrew: return the correlation scores of the row identified by index to all other rows, do not include to itself
		score_distribution = np.delete(corrcoefMat[index],index)

	#print score_distribution

	#if nans != None:
	#	nan_list = np.isnan(score_distribution)
	#	#print nan_list
	#	try:
	#		score_distribution[nan_list] = nans
	#	except TypeError:
	#		pass

	return score_distribution


#kdrew: calculates t-values for array of correlation coefficents
#kdrew: ar is an array of correlation coefficents, n is the number of samples
def tvalue_correlation(ar, n):

	#kdrew: t-value formula
	def tvalue(r,n):
		try:
			#kdrew: get a divide by zero if r is 1.0 (perfect correlation), set r to be close to one to avoid that
			if m.fabs(r-1.0) < EPSILON: 
				r = CLOSE_TO_ONE
				#print "tvalue_correlation: r == 1.0, changing r value to be CLOSE_TO_ONE" 
			if n <= 2:
				#print "tvalue_correlation: degrees of freedom <= 0, returning tvalue = 0.0"
				return 0.0

			t = r/m.sqrt((1-r**2)/(n-2))

		except ZeroDivisionError:
			#print "ZeroDivisionError: r: %s n %s" % (r, n)
			t = 0.0
		except ValueError:
			#print "ValueError: r: %s n %s" % (r, n)
			t = 0.0


		#print "r: %s n: %s t: %s" % (r, n, t)
		return t
		
	tvalue_vectorized = np.vectorize(tvalue)
	return tvalue_vectorized(ar, n)

#kdrew: adding poisson noise might not be the right thing to do here,
#kdrew: essentially we are adding greater variance to larger values
#kdrew: this penalizes larger values too much
def sample_correlation_distribution(matrix, noise_constant=0.0, normalize=False, index=None, iterations=1000, sample_module=None):

	scores_total = None
	tvalues_total = None
	for i in xrange(iterations):
		noise_mat = nu.add_noise(matrix, noise_constant)
		#print "noise_mat: %s" % (noise_mat,)
		sample_mat = nu.sample_noise(noise_mat, sample_module=sample_module)	
		#print "sample_mat: %s" % (sample_mat,)

		if normalize:
			#kdrew: TODO: here we are only normalizing over the bicluster and not the whole data_matrix, should we be normalizing across the whole? probably
			sample_mat = nu.normalize_over_rows(sample_mat)

		scores = correlation_distribution(sample_mat, index=index)
		#print "scores: %s" % (scores,)
		#kdrew: the length parameter is the number of columns in input matrix
		tvalues = tvalue_correlation( scores, matrix.shape[1] )
		#print "tvalues: %s" % (tvalues,)

		if scores_total == None:
			scores_total = scores
		else:
			scores_total += scores

		if tvalues_total == None:
			tvalues_total = tvalues
		else:
			tvalues_total += tvalues

		#print scores_total
		#print tvalues_total


	return ( 1.0*scores_total/iterations, 1.0*tvalues_total/iterations )
	

if __name__ == "__main__":
	main()


