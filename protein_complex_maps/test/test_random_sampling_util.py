#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import numpy as np


class RandomSamplingTest(unittest.TestCase):

	def setUp(self,):

		np.random.seed(123)
		self.data_matrix = np.matrix( np.arange(40).reshape(4,10))

	def testScoreDistribution(self, ):

		dist_dict = rsu.random_sampling_score_distribution( self.data_matrix, su.multiple_dot, rows=[1,2], columns=[0,3], sample_module = np.random )
            
		#print dist_dict
		for t in ['all','columns','rows']:
			print "t: %s, mean: %s, std: %s" % (t, dist_dict[t].mean(), dist_dict[t].std())

		np.testing.assert_almost_equal( dist_dict['all'].mean(), 173.19075 )
		np.testing.assert_almost_equal( dist_dict['all'].std(), 137.794745643 )
		np.testing.assert_almost_equal( dist_dict['columns'].mean(), 183.97625 )
		np.testing.assert_almost_equal( dist_dict['columns'].std(), 37.5867213712 )
		np.testing.assert_almost_equal( dist_dict['rows'].mean(), 116.9875 )
		np.testing.assert_almost_equal( dist_dict['rows'].std(), 119.01008295 )


if __name__ == "__main__":
	unittest.main()


