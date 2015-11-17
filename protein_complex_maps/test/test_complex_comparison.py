#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.complex_comparison as cc
import numpy as np


class ComplexComparisonTest(unittest.TestCase):

    def setUp(self,):

        self.gold_standard = [['a','b','c'],['d','e','f'],['f','g'],['a','b'],['g','i']]
        self.clusters = [['a','b','c'],['d','e','f','h'],['a','b'],['j','k']]

    def testComplexComparison(self, ):
        ccobj = cc.ComplexComparison(self.gold_standard, self.clusters)
        np.testing.assert_almost_equal( ccobj.sensitivity(), 0.75)
        np.testing.assert_almost_equal( ccobj.ppv(), 8.0/13.0)
        np.testing.assert_almost_equal( ccobj.acc(), 0.6793662204867574)

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

if __name__ == "__main__":
        unittest.main()


