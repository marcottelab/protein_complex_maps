#!/usr/bin/python

# Unit tests for FeatureGenerator class

import unittest 
import protein_complex_maps.feature_generator as fg
import protein_complex_maps.bicluster.bicluster as bc
import numpy as np


class FeatureGeneratorTest(unittest.TestCase):

	def setUp(self,):
		print ""
		print np.array([1, 2, 2, 1, 1, 2, 3, 4, 5, 6, 1, 1, 2, 1, 2, 3, 4, 5, 6, 7, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 1, 2, 4, 5, 6, 7, 8, 9]).reshape(4,10)
		self.feature_obj = fg.FeatureGenerator(bc.Bicluster(rows=[0,1],cols=[0,1,2,3]) ,np.array([1, 2, 2, 1, 1, 2, 3, 4, 5, 6, 1, 1, 2, 1, 2, 3, 4, 5, 6, 7, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 1, 2, 4, 5, 6, 7, 8, 9]).reshape(4,10), seed=123 )

	def testRowCorrelationFeatures(self, ):
		self.feature_obj.correlation_feature_row()
		featmat = self.feature_obj.get_row_feature_matrix()

		np.testing.assert_almost_equal( featmat['corr_mean_bc'][0], 0.57735026919 )
		np.testing.assert_almost_equal( featmat['corr_mean_bc'][1], 0.57735026919 )
		np.testing.assert_almost_equal( featmat['corr_mean_bc'][2], -0.788675134595 )
		np.testing.assert_almost_equal( featmat['corr_mean_bc'][3], -0.455341801261 )

		np.testing.assert_almost_equal( featmat['corr_mean_rand'][0], -0.788675134595 )
		np.testing.assert_almost_equal( featmat['corr_mean_rand'][1], -0.455341801261 )
		np.testing.assert_almost_equal( featmat['corr_mean_rand'][2], 0.57735026919 )
		np.testing.assert_almost_equal( featmat['corr_mean_rand'][3], 0.57735026919 )

		np.testing.assert_almost_equal( featmat['corr_mean_ratio'][0], -0.732050807569 )
		np.testing.assert_almost_equal( featmat['corr_mean_ratio'][1], -1.26794919243 )
		np.testing.assert_almost_equal( featmat['corr_mean_ratio'][2], -1.36602540378 )
		np.testing.assert_almost_equal( featmat['corr_mean_ratio'][3], -0.788675134595 )

	def testColumnCorrelationFeatures(self, ):
		self.feature_obj.correlation_feature_column()
		featmat = self.feature_obj.get_column_feature_matrix()
		print featmat
		np.testing.assert_almost_equal( featmat['corr_gain_bc'][0], 0.0773502691896)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][0], -0.0116005290132)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][1], -0.42264973081)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][1], -0.0555270169962)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][2], 0.57735026919)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][2], -0.0165784279392)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][3], 0.0773502691896)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][3], -0.0116005290132)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][4], 0.391124440037)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][4], 0.0)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][5], 0.402787636739)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][5], 0.0)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][6], -0.410683602523)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][6], 0.0)

		np.testing.assert_almost_equal( featmat['corr_gain_bc'][7], 0.0350221665062)
		np.testing.assert_almost_equal( featmat['corr_gain_rand'][7], 0.0)
if __name__ == "__main__":
	unittest.main()


