

import numpy as np

#kdrew: create distribution of scores (score_func) of randomly sampled matrices, 
#kdrew: types is what to randomly sample (rows = rows will be sampled and columns will be fixed)
#kdrew: returns dictionary of distributions keyed by type
def random_sampling_score_distribution( data_matrix, score_func, rows=[], columns=[], types=["all", "rows", "columns"], iterations=1000, sample_module=np.random):

	total_rows, total_columns = data_matrix.shape

	rows_columns_sampled_scores = []
	columns_sampled_scores = []
	rows_sampled_scores = []

	return_dict = dict()

	for t in types:
		return_dict[t] = np.array([])
		for i in xrange(0,iterations):
			if "all" == t:
				#kdrew: random sample from data_matrix keeping size of rows and columns same size as parameters
				sample_rows = sample_module.choice(xrange(0,total_rows), len(rows), replace=False)
				sample_columns = sample_module.choice(xrange(0,total_columns), len(columns), replace=False)

			if "columns" == t:
				#kdrew: random  sample from data_matrix keeping rows the same and columns same size
				sample_rows = rows
				sample_columns = sample_module.choice(xrange(0,total_columns), len(columns), replace=False)

			if "rows" == t:
				#kdrew: random sample from data_matrix keeping columns the same and rows same size 
				sample_rows = sample_module.choice(xrange(0,total_rows), len(rows), replace=False)
				sample_columns = columns

			sample_rows.sort()
			sample_columns.sort()

			submatrix = data_matrix[np.ix_(sample_rows, sample_columns)]
			score = score_func( submatrix ) 
			return_dict[t] = np.append( return_dict[t], score )


	return return_dict
	


