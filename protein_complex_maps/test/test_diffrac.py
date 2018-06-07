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

        d3 = dict()
        d3['ids'] = ['a','b','c','d']
        d3['frac1'] = [1.0,1.0,0.0,1.0]
        d3['frac2'] = [0.0,0.0,0.0,0.0]
        d3['frac3'] = [0.0,0.0,1.0,0.0]
        d3['frac4'] = [0.0,0.0,0.0,0.0]
        df3 = pd.DataFrame(d3, dtype=float)
        df3 = df3.set_index('ids')
        self.elution3 = eff.Elut()
        self.elution3.df = df3

        d4 = dict()
        d4['ids'] = ['e','f','g','h']
        d4['frac1'] = [10.0,10.0,0.0,20.0]
        d4['frac2'] = [0.0,0.0,10.0,40.0]
        d4['frac3'] = [20.0,10.0,0.0,0.0]
        d4['frac4'] = [0.0,0.0,60.0,10.0]
        df4 = pd.DataFrame(d4, dtype=float)
        df4 = df4.set_index('ids')
        self.elution4 = eff.Elut()
        self.elution4.df = df4

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

    def testZscore(self, ):
        print "zscore"
        df = pd.DataFrame()
        mean_abundance = dfc.calc_mean_abundance(self.elution, self.elution2)
        mean_abundance.name = 'mean_abundance'
        df = df.join(mean_abundance, how='outer')

        feature_series = dfc.calc_diffrac(self.elution, self.elution2, normalize_totalCounts=True)
        feature_series.name = 'diffrac_normalized'
        df = df.join(feature_series, how='outer')

        series = pd.Series([True,True,False,False,False],index= ['a','b','c','d','e'])
        series.name = 'annotated'
        df = df.join(series, how='outer')

        print df
        zscores = dfc.calc_zscore(df)

        print zscores.loc['e'] 
        np.testing.assert_almost_equal(zscores.loc['e'],1.4140145096382672)

    def testSlidingZscore(self, ):
        print "sliding zscore"
        df = pd.DataFrame()
        mean_abundance = dfc.calc_mean_abundance(self.elution, self.elution2)
        mean_abundance.name = 'mean_abundance'
        df = df.join(mean_abundance, how='outer')

        feature_series = dfc.calc_diffrac(self.elution, self.elution2, normalize_totalCounts=True)
        feature_series.name = 'diffrac_normalized'
        df = df.join(feature_series, how='outer')

        series = pd.Series([True,True,False,False,False],index= ['a','b','c','d','e'])
        series.name = 'annotated'
        df = df.join(series, how='outer')

        print df
        sliding_zscores = dfc.calc_sliding_zscore(df, window=1)
        print sliding_zscores
        np.testing.assert_almost_equal(sliding_zscores.loc['e']['sliding_zscore'],103.22222222222224)

        sliding_zscores = dfc.calc_sliding_zscore(df, window=1, use_gmm=True)
        print sliding_zscores
        np.testing.assert_almost_equal(sliding_zscores.loc['e']['sliding_zscore'],103.22222222222224)

        print "testing window = 5"
        sliding_zscores = dfc.calc_sliding_zscore(df, window=5)
        print sliding_zscores
        np.testing.assert_almost_equal(sliding_zscores.loc['e']['sliding_zscore'],103.22222222222224)

        print "testing gmm with artifically low min weight threshold 0.2"
        sliding_zscores = dfc.calc_sliding_zscore(df, window=5, use_gmm=True, min_weight_threshold=0.2)
        print sliding_zscores
        #np.testing.assert_almost_equal(sliding_zscores.loc['e']['sliding_zscore'],529.68150807816676) #kdrew: had some instability because both components had equal weight
        np.testing.assert_almost_equal(sliding_zscores.loc['b']['sliding_zscore'],-1.2864445660676826)

    def testSlidingZscore2(self, ):
        print "sliding zscore2"
        df = pd.DataFrame()
        mean_abundance = dfc.calc_mean_abundance(self.elution3, self.elution4)
        mean_abundance.name = 'mean_abundance'
        df = df.join(mean_abundance, how='outer')

        feature_series = dfc.calc_diffrac(self.elution3, self.elution4, normalize_totalCounts=True)
        feature_series.name = 'diffrac_normalized'
        df = df.join(feature_series, how='outer')

        #series = pd.Series([True,True,False,False,False],index= ['a','b','c','d','e'])
        #series.name = 'annotated'
        #df = df.join(series, how='outer')

        print df
        sliding_zscores = dfc.calc_sliding_zscore(df, window=1)

        print sliding_zscores

        #kdrew: a is nan because there is no defined standard deviation within in the window, (b and c are both 1.0
        np.testing.assert_almost_equal(sliding_zscores.loc['a']['sliding_zscore'],np.nan)
        

    def testSlidingZscoreFDR(self, ):
        print "sliding zscore FDR"
        df = pd.DataFrame()
        mean_abundance = dfc.calc_mean_abundance(self.elution, self.elution2)
        mean_abundance.name = 'mean_abundance'
        df = df.join(mean_abundance, how='outer')

        feature_series = dfc.calc_diffrac(self.elution, self.elution2, normalize_totalCounts=True)
        feature_series.name = 'diffrac_normalized'
        df = df.join(feature_series, how='outer')

        series = pd.Series([True,True,False,False,False],index= ['a','b','c','d','e'])
        series.name = 'annotated'
        df = df.join(series, how='outer')

        sliding_zscores = dfc.calc_sliding_zscore(df, window=1)
        df = df.join(sliding_zscores, how='outer')
        fdr_pvals = dfc.calc_sliding_fdr_correct(df)
        print fdr_pvals
        np.testing.assert_almost_equal(fdr_pvals.loc['a']['sliding_pvalues'],0.848460)
        np.testing.assert_almost_equal(fdr_pvals.loc['e']['sliding_pvalues'],0.0)
        np.testing.assert_almost_equal(fdr_pvals.loc['a']['sliding_pvalues_fdrcor'],0.850628)
        np.testing.assert_almost_equal(fdr_pvals.loc['e']['sliding_pvalues_fdrcor'],0.0)

if __name__ == "__main__":
        unittest.main()


