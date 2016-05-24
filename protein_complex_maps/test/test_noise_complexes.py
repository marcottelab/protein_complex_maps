#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.complexes.noise_complexes as nc
import numpy as np
import random
import numpy.random as nr

from scipy.misc import comb


class NoiseComplexesTest(unittest.TestCase):

    def setUp(self,):

        random.seed(1234)
        nr.seed(1234)

        self.clusters = [['a','b','c','d','e','f','g'], 
                ['h','i','j','k','l','m'],
                ['n','o','p','q','r','s','t','u'], 
                ['v','w','x','y','z','aa','ab','ac','ad','ae','af','ag','ao','ap','aq','ar'],
                ['v','w','x','aq','ar','as','at','au']]

        self.clusters2 = [['a','b','c'], 
                ['h','i','j']]

        print "setUp"
        print self.clusters


    def testShuffleComplexes(self, ):
        shuffled_complexes = nc.shuffle_complexes(self.clusters, shuffle_fraction=1.0)
        print "ShuffleComplexes"
        print shuffled_complexes
        #np.testing.assert_almost_equal( ccobj.sensitivity(), 0.75)
        #np.testing.assert_almost_equal( ccobj.ppv(), 8.0/13.0)
        #np.testing.assert_almost_equal( ccobj.acc(), 0.6793662204867574)

    def testBreakupComplexes(self, ):
        broken_complexes = nc.breakup_complexes(self.clusters, breakup_fraction=1.0)
        print "BreakupComplexes"
        print broken_complexes
        assert len(broken_complexes) == 10

    def testBreakupComplexesKeepOriginal(self, ):
        broken_complexes = nc.breakup_complexes(self.clusters, breakup_fraction=1.0, keep_original=True)
        print "BreakupComplexesKeepOriginal"
        print broken_complexes
        assert len(broken_complexes) == 15

    def testRemoveByThreshold(self, ):
        #kdrew: test default, should return same
        remove_by_size_complexes = nc.remove_by_size_complexes(self.clusters)
        print "RemoveByThresholdTest"
        print remove_by_size_complexes
        assert len(remove_by_size_complexes) == len(self.clusters)

    def testRemoveByThreshold2(self, ):
        #kdrew: test lower bound, remove clusters size 7 and below
        remove_by_size_complexes = nc.remove_by_size_complexes(self.clusters, lower_size_threshold=7)
        print "RemoveByThresholdTest2"
        print remove_by_size_complexes
        assert len(remove_by_size_complexes) == 3

    def testPowerset(self, ):
        powerset_comp = nc.powerset_complexes(self.clusters2)
        print "testPowerset"
        ids = set()
        for c in self.clusters2:
            ids = ids.union(set(c))
        #kdrew: calculate total number of combinations that should be in powerset (k>2)
        total_combinations = sum([comb(len(ids),x) for x in range(2,len(ids)+1)])

        print total_combinations
        print len(powerset_comp)
        assert len(powerset_comp) == int(total_combinations)

    def testAddRandomComplex(self, ):
        random_complexes = nc.add_random_complexes(self.clusters, fraction=1.0)
        print "testAddRandomComplex"
        print random_complexes
        assert len(random_complexes) == 10

    def testRemoveFraction(self, ):
        random_complexes = nc.remove_fraction_complexes(self.clusters, fraction=0.5)
        print "testRemoveFraction"
        print random_complexes
        assert len(random_complexes) == 2
        assert len(random_complexes[0]) == 6
        assert len(random_complexes[1]) == 7

    def testRemoveFraction2(self, ):
        random_complexes = nc.remove_fraction_complexes(self.clusters, fraction=0.5, remove_upper=False)
        print "testRemoveFraction2"
        print random_complexes
        assert len(random_complexes) == 3
        assert len(random_complexes[0]) == 8
        assert len(random_complexes[1]) == 8
        assert len(random_complexes[2]) == 16

    def testRemoveSubunitsByComplex(self, ):
        random_complexes,n = nc.remove_subunits_by_complex(self.clusters, fraction=1.0 )
        print "testRemoveSubunitsByComplex"
        print random_complexes
        random_complexes = sorted(random_complexes,key=len)
        assert len(random_complexes) == 5
        assert len(random_complexes[0]) == 5
        assert len(random_complexes[1]) == 6
        assert len(random_complexes[2]) == 7
        assert len(random_complexes[3]) == 7
        assert len(random_complexes[4]) == 15

    def testRemoveSubunits(self, ):
        complexes = sorted(self.clusters,key=len)
        random_complexes = nc.remove_subunits(complexes, totalNumberOfSubunits = 3, shuffle=False )
        print "testRemoveSubunits"
        print random_complexes
        assert len(random_complexes) == 5
        assert len(random_complexes[0]) == 3
        assert len(random_complexes[1]) == 7
        assert len(random_complexes[2]) == 8
        assert len(random_complexes[3]) == 8
        assert len(random_complexes[4]) == 16

    def testRemoveSubunits2(self, ):
        random_complexes = nc.remove_subunits(self.clusters, totalNumberOfSubunits = 10)
        print sum([len(c) for c in random_complexes]) 
        print sum([len(c) for c in self.clusters]) - 10
        assert sum([len(c) for c in random_complexes]) == sum([len(c) for c in self.clusters]) - 10

if __name__ == "__main__":
        unittest.main()


