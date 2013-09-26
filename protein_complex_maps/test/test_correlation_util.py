#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.correlation_util as cu
import numpy as np


class CorrelationTest(unittest.TestCase):

	def setUp(self,):

		self.data_matrix = np.arange(40).reshape(4,10)

		self.array1 = np.array([1,9,2,8,3,7,4,6,5,0])
		self.data_matrix = np.vstack([self.data_matrix, self.array1])

		#print self.data_matrix

	def testMatrix(self, ):
		scores = cu.correlation_distribution(self.data_matrix)
		assert( scores[0] == 1.0 )
		assert( len(scores) == 10 )


	def testArrayVMatrix(self,):
		scores = cu.correlation_distribution(self.data_matrix, 4)
		assert( scores[0] == -0.15151515151515152 )
		assert( len(scores) == 4 )
		


if __name__ == "__main__":
	unittest.main()


