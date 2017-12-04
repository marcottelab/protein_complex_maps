#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.features.ExtractFeatures.Features as eff
import numpy as np
import pandas as pd


class SharedBaitFeatureTest(unittest.TestCase):

    def setUp(self,):

        d = dict()
        d['ids'] = ['a','b','c','d']
        d['frac1'] = [1.0,2.0,0.0,0.0]
        d['frac2'] = [0.0,0.0,3.0,4.0]
        d['frac3'] = [2.0,1.0,0.0,0.0]
        d['frac4'] = [0.0,0.0,6.0,8.0]
        df = pd.DataFrame(d, dtype=float)
        df = df.set_index('ids')
        print df
        print df.dtypes
        self.elution = eff.ElutFeatures(df)
        #feature_matrix = elution.extract_features(feature=args.feature,resampling=args.resampling,iterations=args.iterations,threshold=args.threshold)

    def testSumDifference(self, ):
        sdifs = self.elution.extract_features("sum_difference")
        print sdifs
        assert(sdifs.query("(ID1 == 'a') and (ID2 == 'b')")['sum_difference'].values[0] == 2.0)
        assert(sdifs.query("(ID1 == 'a') and (ID2 == 'c')")['sum_difference'].values[0] == 12.0)
        assert(sdifs.query("(ID1 == 'a') and (ID2 == 'd')")['sum_difference'].values[0] == 15.0)
        assert(sdifs.query("(ID1 == 'c') and (ID2 == 'd')")['sum_difference'].values[0] == 3.0)
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values, 2.302585 )

    def testPearson(self, ):
        pearson_coeffs = self.elution.extract_features("pearsonR")
        print "pearson_coeffs"
        print pearson_coeffs
        #kdrew: hardcoded values are from scipy.stats.pearsonr
        np.testing.assert_almost_equal(pearson_coeffs.query("(ID1 == 'a') and (ID2 == 'b')")['pearsonR'].values[0], 0.63636363636363635)
        np.testing.assert_almost_equal(pearson_coeffs.query("(ID1 == 'a') and (ID2 == 'c')")['pearsonR'].values[0], -0.81818181818181823)
        np.testing.assert_almost_equal(pearson_coeffs.query("(ID1 == 'a') and (ID2 == 'd')")['pearsonR'].values[0], -0.81818181818181823)
        np.testing.assert_almost_equal(pearson_coeffs.query("(ID1 == 'c') and (ID2 == 'd')")['pearsonR'].values[0], 1.0)



if __name__ == "__main__":
        unittest.main()


