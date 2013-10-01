
import numpy as np

class Bicluster(object):

	def __init__(self, rows=[], cols=[]):
		#kdrew: define as sets so we don't have duplicates, there should only be one row index and one column index for each actual row or column
		self.__rows = set(rows)
		self.__cols = set(cols)

	def add_row(self, row_index):
		self.__rows.add(row_index)

	def remove_row(self, row_index):
		self.__rows.remove(row_index)
	
	def add_column(self, col_index):
		self.__cols.add(col_index)

	def remove_column(self, col_index):
		self.__cols.remove(col_index)

	def rows(self, ):
		return list(self.__rows)

	def columns(self, ):
		return list(self.__cols)

	#kdrew: returns submatrix defined by the bicluster
	#kdrew: without_indices creates the submatrix without the row specified by given indices
	def get_submatrix(self, matrix, without_rows=[], without_cols=[]):
		#print self.rows
		#print self.cols
		return matrix[np.ix_(list(self.__rows-set(without_rows)),list(self.__cols-set(without_cols)))]

	#kdrew: get single row with columns defined in bicluster
	def get_row(self, matrix, row, without_cols=[]):
		return matrix[np.ix_(list([row,]),list(self.__cols-set(without_cols)))]

	#kdrew: get single column with rows defined in bicluster
	def get_column(self, matrix, column, without_rows=[]):
		return matrix[np.ix_(list(self.__rows-set(without_rows)),list([column,]))]

	#kdrew: returns rows not in bicluster
	def get_outside_rows(self, matrix):
		all_rows = set(np.arange(matrix.shape[0]))
		#print "all_rows", all_rows
		return all_rows - self.__rows

	#kdrew: returns columns not in bicluster
	def get_outside_cols(self, matrix):
		all_cols = set(np.arange(matrix.shape[1]))
		return all_cols - self.__cols


	def get_random_outside_rows(self, matrix, seed=None, without_rows=[]):
		import random as r
		if seed != None:
			r.seed(seed)
		outside_rows = self.get_outside_rows(matrix) - set(without_rows)
		#print outside_rows
		#kdrew: return a list of random rows the size of the bicluster rows
		random_rows = r.sample(outside_rows,len(self.__rows))
		#print random_rows
		return random_rows

	def get_random_outside_cols(self, matrix, seed=None, without_cols=[]):
		import random as r
		if seed != None:
			r.seed(seed)
		outside_cols = self.get_outside_cols(matrix) - set(without_cols)
		#print outside_cols
		#kdrew: return a list of random cols the size of the bicluster cols
		random_cols = r.sample(outside_cols,len(self.__cols))
		#print random_cols
		return random_cols

