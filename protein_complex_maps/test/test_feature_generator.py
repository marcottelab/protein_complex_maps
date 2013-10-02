#!/usr/bin/python

# Unit tests for FeatureGenerator class

import unittest 
import protein_complex_maps.feature_generator as fg
import protein_complex_maps.bicluster.bicluster as bc
import numpy as np


class FeatureGeneratorTest(unittest.TestCase):

	def setUp(self,):
		np.random.seed(120)
		print ""
		print np.array([1, 2, 2, 1, 1, 2, 3, 4, 5, 6, 1, 1, 2, 1, 2, 3, 4, 5, 6, 7, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 1, 2, 4, 5, 6, 7, 8, 9]).reshape(4,10)
		self.feature_obj = fg.FeatureGenerator(bc.Bicluster(rows=[0,1],cols=[0,1,2,3]) ,np.array([1, 2, 2, 1, 1, 2, 3, 4, 5, 6, 1, 1, 2, 1, 2, 3, 4, 5, 6, 7, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 1, 2, 4, 5, 6, 7, 8, 9]).reshape(4,10), seed=123 )

		self.logreg_feature_obj = fg.FeatureGenerator(bc.Bicluster(rows=[0,1,7],cols=[0,1,2,3]) ,np.array(np.random.random_integers(1,10,80)).reshape(8,10), seed=123 )
		#print np.array(np.random.random_integers(1,10,80)).reshape(8,10)

	def testRowCorrelationFeatures(self, ):
		self.feature_obj.correlation_feature_row()
		featmat = self.feature_obj.get_row_feature_matrix()
		print featmat

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

	def testLogRegression(self, ):
		import pandas as pd
		self.logreg_feature_obj.correlation_feature_row()
		featmat = self.logreg_feature_obj.get_row_feature_matrix()
		print featmat
		test_columns = ['label','corr_mean_bc','corr_mean_rand']
		row_logreg, row_logit_result = self.logreg_feature_obj.create_logistic_regression(featmat, test_columns)
		print row_logit_result.summary()
		print featmat
		new_test_point = pd.DataFrame(np.array([1.0,0.0,1.0]).reshape(1,3), columns=['corr_mean_bc','corr_mean_rand','intercept'])
		print new_test_point
		predictions = row_logit_result.predict(new_test_point)
		print predictions[0]
		#kdrew: CAUTION: this was not done by hand, just assumed it was correct from the run
		np.testing.assert_almost_equal( predictions[0], 0.41736007 )


		bicluster1 = bc.Bicluster(rows=[0,1,7,5],cols=[0,1,2,3])
		bc_corr_mean = self.logreg_feature_obj.correlation_feature_row_by_index(bicluster1, 5)
 		rand_bicluster = self.logreg_feature_obj.get_rand_row_bicluster()
		rand_bicluster.add_row(5)
		rand_bc_corr_mean = self.logreg_feature_obj.correlation_feature_row_by_index(rand_bicluster, 5)
		rand_bicluster.remove_row(5)

		print bc_corr_mean, rand_bc_corr_mean
		new_test_point2 = pd.DataFrame(np.array([bc_corr_mean, rand_bc_corr_mean, 1.0]).reshape(1,3), columns=['corr_mean_bc','corr_mean_rand','intercept'])
		predictions = row_logit_result.predict(new_test_point2)
		print predictions[0]
		#kdrew: CAUTION: this was not done by hand, just assumed it was correct from the run
		np.testing.assert_almost_equal( predictions[0], 0.45052775747620843)

if __name__ == "__main__":
	unittest.main()


