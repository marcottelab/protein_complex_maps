#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.complex_comparison as cc
import numpy as np
import random
import numpy.random as nr


class ComplexComparisonTest(unittest.TestCase):

    def setUp(self,):

        random.seed(1234)
        nr.seed(1234)

        self.gold_standard = [['a','b','c'],['d','e','f'],['f','g'],['a','b'],['g','i']]
        self.clusters = [['a','b','c'],['d','e','f','h'],['a','b'],['j','k']]

        self.gold_standard2 = [['a','b','c','w'],['d','e','f'],['f','g'],['a','b'],['g','i']]
        self.clusters2 = [['a','b','c','w','x','y','z'],['d','e','f','h'],['a','b'],['j','k']]

        self.gold_standard3 = [['a','b','c'],['d','e','f'],['f','g'],['a','b'],['g','i']]
        self.clusters3 = [['a','b','c'],['d','e','f','h'],['a','b'],['j','k'],['c','d'],['f','g','a']]

    def testComplexComparison(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard, self.clusters)
        np.testing.assert_almost_equal( ccobj.sensitivity(), 0.75)
        np.testing.assert_almost_equal( ccobj.ppv(), 8.0/13.0)
        np.testing.assert_almost_equal( ccobj.acc(), 0.6793662204867574)

        print ccobj.sensitivity()
        print ccobj.ppv()
        print ccobj.acc()

    def testComplexComparison2(self, ):
        ccobj2 = cc.ComplexComparison(self.gold_standard2, self.clusters2)
        print ccobj2.sensitivity()
        print ccobj2.ppv()
        print ccobj2.acc()

    def testComplexComparison3(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard, self.clusters)
        df = ccobj.get_na_table()
        np.testing.assert_almost_equal(df[1][1], 9.0/12)
        np.testing.assert_almost_equal(df[0][3], 4.0/6)
        np.testing.assert_almost_equal(ccobj.max_matching_ratio(), 0.575)

    def testComplexComparisonMMR(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard, self.clusters)
        #MMR = sum( maxNA ) / |gs|
        # maxNA(gs[0]) = 3**2 / 3*3 = 1.0
        # maxNA(gs[1]) = 3**2 / 3*4 = 0.75
        # maxNA(gs[2]) = 1**2 / 4*2 = 0.125
        # maxNA(gs[3]) = 2**2 / 2*2 = 1.0
        # maxNA(gs[4]) = 0**2 / 2*2 = 0.0
        # MMR = 2.875 / 5 = 0.575
        print ccobj.mmr()
        np.testing.assert_almost_equal( ccobj.mmr(), 0.575)
        self.clusters.append(['h','k'])
        self.clusters.append(['j','l'])
        ccobj = cc.ComplexComparison(self.gold_standard, self.clusters)
        print ccobj.mmr()
        np.testing.assert_almost_equal( ccobj.mmr(), 0.575)

    def testCliqueComparison(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3)

        result_dict = ccobj.clique_comparison(clique_size=3)
        assert result_dict['tp'] == 6696
        assert result_dict['fp'] == 3304
        #kdrew: actual
        #assert result_dict['tp'] == 2
        #assert result_dict['fp'] == 1
        assert result_dict['fn'] == 0

        result_dict = ccobj.clique_comparison(clique_size=2)
        assert result_dict['tp'] == 7309
        assert result_dict['fp'] == 2691
        assert result_dict['fn'] == 1059
        #kdrew: actual
        #assert result_dict['tp'] == 8
        #assert result_dict['fp'] == 3
        #assert result_dict['fn'] == 1

    def testCliqueComparisonMetric(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3)
        d = ccobj.clique_comparison_metric()

        #kdrew: actual
        np.testing.assert_almost_equal( d[2]['precision'], 0.7272727272727273, 2)
        #kdrew: estimated
        np.testing.assert_almost_equal( d[2]['precision'], 0.7272)

        #kdrew: actual
        np.testing.assert_almost_equal( d[2]['recall'], 0.8888888888888888, 1 )
        #kdrew: estimated
        #kdrew: old value before bug fix on recall calculation
        #np.testing.assert_almost_equal( d[2]['recall'], 0.8678839957035446)
        np.testing.assert_almost_equal( d[2]['recall'], 0.8893)

        #kdrew: actual
        np.testing.assert_almost_equal( d[3]['precision'], 0.666666666667, 1 )
        #kdrew: estimated
        np.testing.assert_almost_equal( d[3]['precision'], 0.6771 )

        #kdrew: actual
        np.testing.assert_almost_equal( d[3]['recall'], 1.0 )
        #kdrew: estimated
        np.testing.assert_almost_equal( d[3]['recall'], 1.0 )

        ccmm = ccobj.clique_comparison_metric_mean()
        #kdrew: old value before bug fix on recall calculation
        #np.testing.assert_almost_equal(ccmm['recall_mean'], 0.93259774526265293)
        #kdrew: (0.8876 + 1.0) / 2
        np.testing.assert_almost_equal(ccmm['recall_mean'], 0.94379999999999997)
        np.testing.assert_almost_equal(ccmm['precision_mean'], 0.69205000000000005)



if __name__ == "__main__":
        unittest.main()


