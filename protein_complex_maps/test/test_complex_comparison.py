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

        self.gold_standard4 = [['a','b','c','d','e','f'],['g','h','i','j','k','l','m','n'],['o','p'],['q','r'],['s','t'],['u','v']]
        #self.clusters4 = [['a','b','c','d','e','f'],['g','h','i','j','k','l','m','n'],['a','g','o','q','s','u']]
        self.clusters4 = [['a','b','c','d','e','f'],['g','h','i','j','k','l'],['a','g','o','q','s','u']]

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
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0)

        result_dict = ccobj.clique_comparison(clique_size=3)
        print "testCliqueComparison: %s" % (result_dict,)
        #kdrew: previous way of sampling, new way below removes redudancy
        #np.testing.assert_equal(result_dict['tp'],6696)
        #np.testing.assert_equal(result_dict['fp'],3304)
        np.testing.assert_equal(result_dict['tp'],2)
        np.testing.assert_equal(result_dict['fp'],1)
        #kdrew: actual
        #assert result_dict['tp'] == 2
        #assert result_dict['fp'] == 1
        np.testing.assert_equal(result_dict['fn'],0)

        #kdrew: [a,b,c] and [d,e,f] are both in gold standard so both are true positives, 
        #kdrew: [d,e,h],[e,f,h] and [d,f,h] are all ignored because 'h' does not exist in gold standard at all
        #kdrew: [f,g,a] is not in gold standard but all of the proteins are so it is a false positive


        result_dict = ccobj.clique_comparison(clique_size=2)
        #kdrew: previous way of sampling, new way below removes redudancy
        #np.testing.assert_equal (result_dict['tp'], 7309)
        #np.testing.assert_equal (result_dict['fp'], 2691)
        #np.testing.assert_equal (result_dict['fn'], 1059)
        np.testing.assert_equal (result_dict['tp'], 7)
        np.testing.assert_equal (result_dict['fp'], 3)
        np.testing.assert_equal (result_dict['fn'], 1)
        #kdrew: actual
        #assert result_dict['tp'] == 8 #kdrew: tp == 8 if we *are not* removing redundant cliques (ab is listed twice)
        #assert result_dict['tp'] == 7 #kdrew: tp == 7 if we *are* removing  redundant cliques
        #assert result_dict['fp'] == 3
        #assert result_dict['fn'] == 1

    def testCliqueComparisonMetricExact(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, exact=True)
        d_exact = ccobj.clique_comparison_metric()
        np.testing.assert_almost_equal( d_exact[2]['precision'], 0.7272727272727273)
        np.testing.assert_almost_equal( d_exact[2]['recall'], 0.8888888888888888)
        np.testing.assert_almost_equal( d_exact[3]['precision'], 0.666666666667)
        np.testing.assert_almost_equal( d_exact[3]['recall'], 1.0 )

    def testCliqueComparisonMetric(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0)
        d = ccobj.clique_comparison_metric()

        #kdrew: changed code to remove redundancy 8.0/11.0
        #np.testing.assert_almost_equal( d[2]['precision'], 0.7272727272727273, 2)
        #np.testing.assert_almost_equal( d[2]['precision'], 0.7272)
        #kdrew: actual 7.0/10.0
        np.testing.assert_almost_equal( d[2]['precision'], 0.7, 2)
        #kdrew: estimated
        np.testing.assert_almost_equal( d[2]['precision'], 0.7)

        #kdrew: actual
        np.testing.assert_almost_equal( d[2]['recall'], 0.875, 1 )
        #kdrew: estimated
        #kdrew: old value before bug fix on recall calculation
        #np.testing.assert_almost_equal( d[2]['recall'], 0.8678839957035446)
        #kdrew: before removing redudancy
        #np.testing.assert_almost_equal( d[2]['recall'], 0.8893)
        np.testing.assert_almost_equal( d[2]['recall'], 0.875)

        #kdrew: actual
        np.testing.assert_almost_equal( d[3]['precision'], 0.666666666667, 1 )
        #kdrew: estimated
        #kdrew: old value before removing redundancy
        #np.testing.assert_almost_equal( d[3]['precision'], 0.6771 )
        np.testing.assert_almost_equal( d[3]['precision'], 0.666666666667)

        #kdrew: actual
        np.testing.assert_almost_equal( d[3]['recall'], 1.0 )
        #kdrew: estimated
        np.testing.assert_almost_equal( d[3]['recall'], 1.0 )

        print "f1score 2: %s" % d[2]['f1score']
        print "f1score 3: %s" % d[3]['f1score']
        #np.testing.assert_almost_equal( d[3]['f1score'], 0.80746526742591385 )
        #np.testing.assert_almost_equal( d[2]['f1score'], 0.80012243736467659 )

        np.testing.assert_almost_equal( d[3]['f1score'], 0.80000000000024007)
        np.testing.assert_almost_equal( d[2]['f1score'], 0.7777777777777779)

    

        ccmm = ccobj.clique_comparison_metric_mean()
        #kdrew: old value before bug fix on recall calculation
        #np.testing.assert_almost_equal(ccmm['recall_mean'], 0.93259774526265293)
        #kdrew: (0.8876 + 1.0) / 2
        #np.testing.assert_allclose(ccmm['recall_mean'], 0.94379999999999997, 0.1)
        #np.testing.assert_allclose(ccmm['precision_mean'], 0.69205000000000005, 0.1)


        #kdrew: old values before reduction of redundancy
        #np.testing.assert_allclose(ccmm['precision_mean'], 0.70215)
        #np.testing.assert_allclose(ccmm['recall_mean'], 0.94465)
        #kdrew: (0.7 + 0.666666666667) / 2
        np.testing.assert_allclose(ccmm['precision_mean'], 0.6833333333335)
        #kdrew: (0.875 + 1.0) / 2
        np.testing.assert_allclose(ccmm['recall_mean'], 0.9375)


        grandf1score = ccobj.clique_comparison_metric_grandf1score()
        #kdrew: old value before reduction in redundancy
        #np.testing.assert_allclose(grandf1score, 0.80377708281154325)
        np.testing.assert_allclose(grandf1score, 0.78873239436631382)

    def testCliqueComparisonMetricSampling(self, ):
        print "Sampling test"
        actual_precision = 8.0/11
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=10)
        d = ccobj.clique_comparison_metric()
        print "10 precision: %s" % d[2]['precision']
        print "10 diff: %s" % abs(d[2]['precision'] - actual_precision)
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=100)
        d = ccobj.clique_comparison_metric()
        print "100 precision: %s" % d[2]['precision']
        print "100 diff: %s" % abs(d[2]['precision'] - actual_precision)
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=1000)
        d = ccobj.clique_comparison_metric()
        print "1k precision: %s" % d[2]['precision']
        print "1k diff: %s" % abs(d[2]['precision'] - actual_precision)
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=10000)
        d = ccobj.clique_comparison_metric()
        print "10k precision: %s" % d[2]['precision']
        print "10k diff: %s" % abs(d[2]['precision'] - actual_precision)
        ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=100000)
        d = ccobj.clique_comparison_metric()
        print "100k precision: %s" % d[2]['precision']
        print "100k diff: %s" % abs(d[2]['precision'] - actual_precision)
        #ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=1000000)
        #d = ccobj.clique_comparison_metric()
        #print "1m precision: %s" % d[2]['precision']
        #print "1m diff: %s" % abs(d[2]['precision'] - actual_precision)
        #ccobj = cc.ComplexComparison(self.gold_standard3, self.clusters3, pseudocount=0, samples=10000000)
        #d = ccobj.clique_comparison_metric()
        #print "10m precision: %s" % d[2]['precision']
        #print "10m diff: %s" % abs(d[2]['precision'] - actual_precision)

    def testCliqueComparisonMetricSampling2(self, ):
        print "Sampling test 2"
        ccobj = cc.ComplexComparison(self.gold_standard4, self.clusters4, pseudocount=0)
        d = ccobj.clique_comparison_metric()
        print "precision 2: %s" % d[2]['precision']
        print "precision 3: %s" % d[3]['precision']
        print "precision 4: %s" % d[4]['precision']
        print "precision 5: %s" % d[5]['precision']
        print "precision 6: %s" % d[6]['precision']
        #print "precision 7: %s" % d[7]['precision']
        #print "precision 8: %s" % d[8]['precision']

        ccobj = cc.ComplexComparison(self.gold_standard4, self.clusters4, pseudocount=0, exact=True)
        d_exact = ccobj.clique_comparison_metric()
        print "precision exact 2: %s" % d_exact[2]['precision']
        print "precision exact 3: %s" % d_exact[3]['precision']
        print "precision exact 4: %s" % d_exact[4]['precision']
        print "precision exact 5: %s" % d_exact[5]['precision']
        print "precision exact 6: %s" % d_exact[6]['precision']
        #print "precision exact 7: %s" % d_exact[7]['precision']
        #print "precision exact 8: %s" % d_exact[8]['precision']

        print "diff 2: %s" % abs(d[2]['precision'] - d_exact[2]['precision'])
        print "diff 3: %s" % abs(d[3]['precision'] - d_exact[3]['precision'])
        print "diff 4: %s" % abs(d[4]['precision'] - d_exact[4]['precision'])
        print "diff 5: %s" % abs(d[5]['precision'] - d_exact[5]['precision'])
        print "diff 6: %s" % abs(d[6]['precision'] - d_exact[6]['precision'])
        #print "diff 7: %s" % abs(d[7]['precision'] - d_exact[7]['precision'])
        #print "diff 8: %s" % abs(d[8]['precision'] - d_exact[8]['precision'])

if __name__ == "__main__":
        unittest.main()


