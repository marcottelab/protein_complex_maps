

import numpy as np


#kdrew: create distribution of scores (score_func) of randomly sampled matrices, 
class RandomSampling(object):
	
	def __init__( self, score_function, iterations=1000, sample_module=np.random ):
		self.score_function = score_function
		self.iterations = iterations
		self.sample_module = sample_module
		self.all_cache = {}

	

	#kdrew: add function to cache mean and std for 'all'
	def score_mean_std_all( self, data_matrix, numrows, numcolumns, force=False ):
		try: 
			if force:
				raise KeyError
			else:
				mean, std = self.all_cache[ numrows, numcolumns ]
		except KeyError:
			dist = self.random_sampling_score_distribution_all( data_matrix, numrows, numcolumns )
			mean = dist.mean()
			std = dist.std()
			self.all_cache[numrows, numcolumns] = (mean, std)

		return mean, std

	#kdrew: convience function for printing zscores
	def print_score( self, data_matrix, bicluster ):
		score = self.score_function( bicluster.get_submatrix( data_matrix ) )

		zscore_all = self.zscore_all( data_matrix, bicluster )
		zscore_cols = self.zscore_columns( data_matrix, bicluster )
		zscore_rows = self.zscore_rows( data_matrix, bicluster )

		print "score %s" % (score, )
		print "zscore_all %s" % (zscore_all, )
		print "zscore_cols %s" % (zscore_cols, )
		print "zscore_rows %s" % (zscore_rows, )
		print "zscore_multiplied %s" % (zscore_all * zscore_cols * zscore_rows)

	def zscore_all( self, data_matrix, bicluster ):
		submatrix = bicluster.get_submatrix(data_matrix)
		mean, std = self.score_mean_std_all( data_matrix, numrows=submatrix.shape[0], numcolumns=submatrix.shape[1] )  
		score = self.score_function( submatrix )
		zscore = abs(mean - score)/std
		return zscore

	def zscore_columns( self, data_matrix, bicluster ):
		submatrix = bicluster.get_submatrix(data_matrix)
		distribution = self.random_sampling_score_distribution_columns( data_matrix, rows=bicluster.rows(), columns=bicluster.columns() )
		mean = distribution.mean()
		std = distribution.std()
		score = self.score_function( submatrix )
		zscore = abs(mean - score)/std
		return zscore

	def zscore_rows( self, data_matrix, bicluster ):
		submatrix = bicluster.get_submatrix(data_matrix)
		distribution = self.random_sampling_score_distribution_rows( data_matrix, rows=bicluster.rows(), columns=bicluster.columns() )
		mean = distribution.mean()
		std = distribution.std()
		score = self.score_function( submatrix )
		zscore = abs(mean - score)/std
		return zscore


	def random_sampling_score_distribution_all(self, data_matrix, numrows, numcolumns ):

		total_rows, total_columns = data_matrix.shape

		return_distribution = np.array([])

		for i in xrange(0,self.iterations):
			#kdrew: random sample from data_matrix keeping size of rows and columns same size as parameters
			sample_rows = self.sample_module.choice(xrange(0,total_rows), numrows, replace=False)
			sample_columns = self.sample_module.choice(xrange(0,total_columns), numcolumns, replace=False)

			sample_rows.sort()
			sample_columns.sort()

			submatrix = data_matrix[np.ix_(sample_rows, sample_columns)]
			score = self.score_function( submatrix ) 
			return_distribution = np.append( return_distribution, score )

		return return_distribution

	def random_sampling_score_distribution_columns(self, data_matrix, rows=[], columns=[] ):

		if len(rows) == 0: 
			raise EmptyRowError

		if len(columns) == 0: 
			raise EmptyColumnsError

		total_rows, total_columns = data_matrix.shape

		return_distribution = np.array([])

		for i in xrange(0,self.iterations):
			#kdrew: random  sample from data_matrix keeping rows the same and columns same size
			sample_rows = rows
			sample_columns = self.sample_module.choice(xrange(0,total_columns), len(columns), replace=False)

			sample_rows.sort()
			sample_columns.sort()

			submatrix = data_matrix[np.ix_(sample_rows, sample_columns)]
			score = self.score_function( submatrix ) 
			return_distribution = np.append( return_distribution, score )

		return return_distribution

	def random_sampling_score_distribution_rows(self, data_matrix, rows=[], columns=[] ):

		if len(rows) == 0: 
			raise EmptyRowError

		if len(columns) == 0: 
			raise EmptyColumnsError

		total_rows, total_columns = data_matrix.shape

		return_distribution = np.array([])

		for i in xrange(0,self.iterations):
			#kdrew: random sample from data_matrix keeping columns the same and rows same size 
			sample_rows = self.sample_module.choice(xrange(0,total_rows), len(rows), replace=False)
			sample_columns = columns

			sample_rows.sort()
			sample_columns.sort()

			submatrix = data_matrix[np.ix_(sample_rows, sample_columns)]
			score = self.score_function( submatrix ) 
			return_distribution = np.append( return_distribution, score )

		return return_distribution


#kdrew: child class which holds data_matrix for scoring
class RandomSamplingScore(RandomSampling):
	
	def __init__( self, data_matrix, score_function, iterations=1000, sample_module=np.random ):
		super(RandomSamplingScore, self).__init__(score_function, iterations=iterations, sample_module=sample_module)
		self.data_matrix = data_matrix

	def zscore_all(self, matrix):
		mean, std = self.score_mean_std_all( self.data_matrix, numrows=matrix.shape[0], numcolumns=matrix.shape[1] )  
		score = self.score_function( matrix )
		zscore = abs(mean - score)/std
		return zscore

	def zscore_all_neg(self, matrix):
		return -1.0*self.zscore_all(matrix)
		
class EmptyRowError(Exception):
	def __init__(self, expr, msg):
		self.expr = expr
		self.msg = msg

class EmptyColumnError(Exception):
	def __init__(self, expr, msg):
		self.expr = expr
		self.msg = msg





