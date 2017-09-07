#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.score_util as su
import protein_complex_maps.normalization_util as nu
import numpy as np


class ScoreTest(unittest.TestCase):

	def setUp(self,):

		self.data_matrix = np.arange(40).reshape(4,10)

		#self.array1 = np.array([1,9,2,8,3,7,4,6,5,0])
		#self.data_matrix = np.vstack([self.data_matrix, self.array1])

		print self.data_matrix

	def testMultipleDot(self, ):
		score = su.multiple_dot_per_unit(self.data_matrix)
		print score
		#assert( score == 72033.3 ) #kdrew: divided by number of columns
		assert( score == 18008.325 ) #kdrew: divided by number of columns and rows

	def testLogMultipleDot(self,):
		dmat = nu.add_noise_over_columns(self.data_matrix)
		#print dmat
		dmat = nu.normalize_over_rows(dmat)
		print dmat
		#log_dmat = nu.log_transform(dmat)
		#print log_dmat
		score = su.sum_cells(dmat)
		print score
		
		#np.testing.assert_almost_equal( score, -1.73220495821 )
		np.testing.assert_almost_equal( score, 0.000343732988528)


if __name__ == "__main__":
	unittest.main()


