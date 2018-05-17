#!/usr/bin/python

# Unit tests for diffrac

import unittest
import numpy as np
import pandas as pd

import scipy.stats as stats

import protein_complex_maps.features.ExtractFeatures.Features as eff
import protein_complex_maps.features.diffrac as dfc


class DifFracTest(unittest.TestCase):

    def setUp(self,):

        d = dict()
        d['ids'] = ['a','b','c','d','e']
        d['frac1'] = [1.0,2.0,0.0,0.0,0.0]
        d['frac2'] = [0.0,0.0,3.0,4.0,0.0]
        d['frac3'] = [2.0,1.0,0.0,0.0,7.0]
        d['frac4'] = [0.0,0.0,6.0,8.0,10.0]
        df = pd.DataFrame(d, dtype=float)
        df = df.set_index('ids')
        self.elution = eff.Elut()
        self.elution.df = df

        d2 = dict()
        d2['ids'] = ['a','b','c','d']
        d2['frac1'] = [1.0,1.0,0.0,2.0]
        d2['frac2'] = [0.0,0.0,1.0,4.0]
        d2['frac3'] = [2.0,1.0,0.0,0.0]
        d2['frac4'] = [0.0,0.0,6.0,10.0]
        df2 = pd.DataFrame(d2, dtype=float)
        df2 = df2.set_index('ids')
        self.elution2 = eff.Elut()
        self.elution2.df = df2

    def testDifFrac(self, ):
        #kdrew: subtracts one dataframe from another and sums the absolute values (L1 norm)
        diffrac_sum = dfc.calc_diffrac(self.elution, self.elution2)
        print diffrac_sum
        assert(diffrac_sum.loc['a'] == 0.0)
        assert(diffrac_sum.loc['b'] == 1.0)
        assert(diffrac_sum.loc['c'] == 2.0)
        assert(diffrac_sum.loc['d'] == 4.0)
        assert(diffrac_sum.loc['e'] == 17.0)

        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values, 2.302585 )

    def testEMD(self, ):
        print "testEMD"
        #kdrew: calculates earth mover's distance
        emd_scores = dfc.calc_emd(self.elution, self.elution2)
        print emd_scores
        assert(emd_scores.loc['a'] == 0.0)
        assert(emd_scores.loc['b'] == 1.0)
        assert(emd_scores.loc['c'] == 2.0)
        assert(emd_scores.loc['d'] == 4.0)
        assert(emd_scores.loc['e'] == 17.0)

        emd_scores = dfc.calc_emd(self.elution2, self.elution)
        print emd_scores

    def testDifCorrelation(self, ):
        correlations = dfc.calc_correlation(self.elution, self.elution2, correlation_func=lambda x,y: stats.pearsonr(x,y)[0])
        print correlations
        assert(correlations.loc['a'] == 1.0)
        np.testing.assert_almost_equal(correlations.loc['b'], 0.904534)
        np.testing.assert_almost_equal(correlations.loc['c'], 0.939394, decimal=5)
        np.testing.assert_almost_equal(correlations.loc['d'], 0.966988, decimal=5)
        assert(correlations.loc['e'] == 0.0)

    def testMeanAbundance(self, ):
        mean_abundance = dfc.calc_mean_abundance(self.elution, self.elution2)
        print mean_abundance
        assert(mean_abundance.loc['a'] == 3.0)
        np.testing.assert_almost_equal(mean_abundance.loc['b'], 2.5)
        np.testing.assert_almost_equal(mean_abundance.loc['c'], 8.0)
        np.testing.assert_almost_equal(mean_abundance.loc['d'], 14.0)
        assert(mean_abundance.loc['e'] == 8.5)

if __name__ == "__main__":
        unittest.main()


