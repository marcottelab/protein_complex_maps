#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.correlation_util as cu
import numpy as np


class CorrelationTest(unittest.TestCase):

	def setUp(self,):

		self.data_matrix = np.arange(50).reshape(5,10)
		#print self.data_matrix

		self.array1 = np.array([1,9,2,8,3,7,4,6,5,0])

	def testMatrix(self, ):
		scores, pvals = cu.correlation_distribution(self.data_matrix)
		assert( scores[0] == 1.0 )
		assert( pvals[0] == 0.0 )
		assert( len(scores) == 10 )


	def testArrayVMatrix(self,):
		scores, pvals = cu.correlation_distribution(self.data_matrix, self.array1)
		assert( scores[0] == -0.15151515151515152 )
		assert( pvals[0] == 0.67606517599785354 )
		assert( len(scores) == 5 )
		


if __name__ == "__main__":
	unittest.main()


