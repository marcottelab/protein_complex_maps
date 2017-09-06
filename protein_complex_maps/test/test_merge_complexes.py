#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.complex_merge as cm
import numpy as np
import random
import numpy.random as nr

from scipy.misc import comb


class MergeComplexesTest(unittest.TestCase):

    def setUp(self,):

        random.seed(1234)
        nr.seed(1234)

        self.clusters = [frozenset(['a','b','c','d','e']), frozenset(['a','b','c']), frozenset(['a','b']), frozenset(['f','g','h','i']), frozenset(['f','g','h','i','j']), frozenset(['k','l','m']), frozenset(['k','l'])]

        print "setUp"
        print self.clusters

        self.complex_size = 4
        self.remove_large_subcomplexes = True
        self.merge_threshold = 0.6
        self.remove_largest = True


    def testShuffleComplexes(self, ):
        complexes = cm.merge_complexes(self.clusters, self)
        print "Merged Complexes"
        print complexes
        assert complexes  == [frozenset(['k','l'])]
        #np.testing.assert_almost_equal( ccobj.sensitivity(), 0.75)
        #np.testing.assert_almost_equal( ccobj.ppv(), 8.0/13.0)
        #np.testing.assert_almost_equal( ccobj.acc(), 0.6793662204867574)


if __name__ == "__main__":
        unittest.main()


