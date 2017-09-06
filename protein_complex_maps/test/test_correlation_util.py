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
		print scores
		assert( scores[0] == 1.0 )
		assert( len(scores) == 10 )


	def testArrayVMatrix(self,):
		scores = cu.correlation_distribution(self.data_matrix, 4)
		assert( scores[0] == -0.15151515151515152 )
		assert( len(scores) == 4 )
		

	def testTvalueCorrelation(self,):
		scores = cu.correlation_distribution(self.data_matrix)
		#print "ttc scores: %s" % (scores,)
		tvalues = cu.tvalue_correlation(scores, 10)
		#print "ttc tvalues: %s" % (tvalues,)
		np.testing.assert_almost_equal( tvalues[0], 19999.9998053 )
		np.testing.assert_almost_equal( tvalues[3], -0.433554984762 )

	def testPoissonCorrelation(self,):

		np.random.seed(12345)
		mat1 = np.array([1,2,3,4,5,6,7,8,9,10]).reshape(2,5)
		mat2 = np.array([10,20,30,40,50,60,70,80,90,100]).reshape(2,5)
		#print mat1
		#print mat2

		scores, tvalues = cu.sample_correlation_distribution(mat1, iterations=1000, sample_module = np.random.poisson)
		#print "tpc scores: %s" % (scores,)
		#print "tpc tvalues: %s" % (tvalues,)
		np.testing.assert_almost_equal( scores[0], 0.30374264 )
		np.testing.assert_almost_equal( tvalues[0], 0.91128502 )

		scores, tvalues = cu.sample_correlation_distribution(mat2, iterations=1000, sample_module = np.random.poisson)
		#print "tpc scores: %s" % (scores,)
		#print "tpc tvalues: %s" % (tvalues,)
		np.testing.assert_almost_equal( scores[0], 0.84185242 )
		np.testing.assert_almost_equal( tvalues[0], 4.0016898 )


	def testNormalCorrelation(self,):
		print "testNormalCorrelation"
		np.random.seed(12345)
		mat1 = np.array([1,2,3,4,5,6,7,8,9,10]).reshape(2,5)
		#mat2 = np.array([10,20,30,40,50,60,70,80,90,100]).reshape(2,5)
		mat2 = np.array([101,102,103,104,105,106,107,108,109,110]).reshape(2,5)
		#print mat1
		#print mat2

		#kdrew: I did not hand calculate these 

		scores, tvalues = cu.sample_correlation_distribution( mat1, iterations=1000, sample_module = np.random.normal )
		#print "tpc scores: %s" % (scores,)
		#print "tpc tvalues: %s" % (tvalues,)
		np.testing.assert_almost_equal( scores[0], 0.73921064 )
		np.testing.assert_almost_equal( tvalues[0], 2.80124664 )

		scores, tvalues = cu.sample_correlation_distribution(mat2, iterations=1000, sample_module = np.random.normal)
		#print "tpc scores: %s" % (scores,)
		#print "tpc tvalues: %s" % (tvalues,)
		np.testing.assert_almost_equal( scores[0], 0.73897173 )
		np.testing.assert_almost_equal( tvalues[0], 2.90779262 )


if __name__ == "__main__":
	unittest.main()


