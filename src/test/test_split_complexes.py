#!/usr/bin/python

# Unit tests for splitting gold standard complexes into test and training

import unittest
import protein_complex_maps.features.split_complexes as sc
import numpy as np
import pandas as pd


class SplitComplexesTest(unittest.TestCase):

    def setUp(self,):

        self.complexes = [['a','b','c'],['c','d','e'],['f','g','h'],['g','h','i'],['a','b'],['d','e'],['b','d'],['a','i']]
        self.complexes1 = [['a','b','c'],['c','d','e'],['f','g','h'],['g','h','i']]
        self.complexes2 = [['a','b'],['d','e'],['b','d'],['a','i']]
        np.random.seed(123)

    def testSplitComplexes(self, ):

        test_complexes, train_complexes, test_ppis, train_ppis, neg_test_ppis, neg_train_ppis = sc.split_complexes(self.complexes)
        print "neg_test_ppis: "
        print neg_test_ppis
        print "neg_train_ppis: "
        print neg_train_ppis

        for ppi in test_ppis:
            assert ppi not in neg_test_ppis
            assert ppi not in train_ppis
            assert ppi not in neg_train_ppis

        for ppi in train_ppis:
            assert ppi not in neg_train_ppis
            assert ppi not in test_ppis
            assert ppi not in neg_test_ppis

        for ppi in neg_train_ppis:
            assert ppi not in train_ppis
            assert ppi not in test_ppis
            assert ppi not in neg_test_ppis

        for ppi in neg_test_ppis:
            assert ppi not in train_ppis
            assert ppi not in test_ppis
            assert ppi not in neg_train_ppis

    def testSplitComplexesRemoveOverlapping(self, ):
        clist1, clist2 = sc.remove_overlapping_complexes(self.complexes1, self.complexes2)
        print "testSplitComplexesRemoveOverlapping"
        print clist1
        print clist2

        #kdrew: check to see if for the overlapping complexes, atleast one of them has been removed
        assert set(['a','b','c']) not in [set(x) for x in clist1] or set(['a','b']) not in [set(x) for x in clist2]
        assert set(['c','d','e']) not in [set(x) for x in clist1] or set(['d','e']) not in [set(x) for x in clist2]
        #kdrew: check to make sure non-overlapping complexes have not been removed
        assert set(['g','h','i']) in [set(x) for x in clist1] 
        assert set(['a','i']) in [set(x) for x in clist2] 

if __name__ == "__main__":
        unittest.main()


