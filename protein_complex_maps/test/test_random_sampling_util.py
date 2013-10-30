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

		mdn_randsamp = rsu.RandomSampling( su.multiple_dot_per_unit, sample_module = np.random )
		dist_all = mdn_randsamp.random_sampling_score_distribution_all( self.data_matrix, numrows=2, numcolumns=2 )
		dist_columns = mdn_randsamp.random_sampling_score_distribution_columns( self.data_matrix, rows=[1,2], columns=[0,3] )
		dist_rows = mdn_randsamp.random_sampling_score_distribution_rows( self.data_matrix, rows=[1,2], columns=[0,3] )

		np.testing.assert_almost_equal( dist_all.mean(), 173.19075 )
		np.testing.assert_almost_equal( dist_all.std(), 137.794745643 )
		np.testing.assert_almost_equal( dist_columns.mean(), 182.19475 )
		np.testing.assert_almost_equal( dist_columns.std(), 37.544273996676246 )
		np.testing.assert_almost_equal( dist_rows.mean(), 112.44 )
		np.testing.assert_almost_equal( dist_rows.std(), 116.115035202 )

	def testAllCache(self, ):
		mdn_randsamp = rsu.RandomSampling( su.multiple_dot_per_unit, sample_module = np.random )
		

		#kdrew: calc for first time
		mean1, std1 = mdn_randsamp.score_mean_std_all( self.data_matrix, 2, 2 )
		#kdrew: grab from cache second time
		mean2, std2 = mdn_randsamp.score_mean_std_all( self.data_matrix, 2, 2 )
		#kdrew: force calc third time
		mean3, std3 = mdn_randsamp.score_mean_std_all( self.data_matrix, 2, 2, force=True )

		#kdrew: grab from cache should be same
		assert(mean1 == mean2)
		assert(std1 == std2)
		np.testing.assert_almost_equal( mean1, mean2 )
		np.testing.assert_almost_equal( std1, std2 )

		#kdrew: force calc should be different
		assert(mean1 != mean3)
		assert(std1 != std3)


if __name__ == "__main__":
	unittest.main()


