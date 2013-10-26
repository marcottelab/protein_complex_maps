

import numpy as np

#kdrew: create distribution of scores (score_func) of randomly sampled matrices, 
#kdrew: types is what to randomly sample (rows = rows will be sampled and columns will be fixed)
#kdrew: returns dictionary of distributions keyed by type

class RandomSampling(object):
	
	def __init__( self, score_function, iterations=1000, sample_module=np.random ):
		self.score_function = score_function
		self.iterations = iterations
		self.sample_module = sample_module
		self.all_cache = {}

	#kdrew: add function to cache mean and std for 'all'
	def score_mean_std_all( self, data_matrix, numrows, numcolumns ):
		try: 
			mean, std = self.all_cache[ numrows, numcolumns ]
		except KeyError:
			dist = self.random_sampling_score_distribution_all( data_matrix, numrows, numcolumns )
			mean = dist.mean()
			std = dist.std()
			self.all_cache[numrows, numcolumns] = (mean, std)

		return mean, std


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


		


