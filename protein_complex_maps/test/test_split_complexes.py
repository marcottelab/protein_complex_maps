#!/usr/bin/python

# Unit tests for splitting gold standard complexes into test and training

import unittest
import protein_complex_maps.preprocessing_util.complexes.split_complexes as sc
import numpy as np
import pandas as pd
import itertools as it


class SplitComplexesTest(unittest.TestCase):

    def setUp(self,):

        self.complexes = [['a','b','c'],['c','d','e'],['f','g','h'],['g','h','i'],['a','b'],['d','e'],['b','d'],['a','i']]
        self.complexes1 = [['a','b','c'],['c','d','e'],['f','g','h'],['g','h','i']]
        self.complexes2 = [['a','b'],['d','e'],['b','d'],['a','i']]
        np.random.seed(123)

    def testSplitComplexes(self, ):

        test_complexes, train_complexes, test_ppis, train_ppis, neg_test_ppis, neg_train_ppis = sc.split_complexes(self.complexes)
        print "test_complexes:"
        print test_complexes
        print "train_complexes:"
        print train_complexes
        print "test_ppis: "
        print test_ppis
        print "train_ppis: "
        print train_ppis
        print "neg_test_ppis: "
        print neg_test_ppis
        print "neg_train_ppis: "
        print neg_train_ppis
        print ""

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

        #kdrew: check to see if there is overlap between test_ppis and training_complexes and vis versa
        for test_ppi in test_ppis:
            assert test_ppi not in [set(x) for train_complex in train_complexes for x in it.combinations(train_complex,2) ] 
        for train_ppi in train_ppis:
            assert train_ppi not in [set(x) for test_complex in test_complexes for x in it.combinations(test_complex,2) ]  


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

    def testSplitComplexes_thresholding(self, ):

        print ""
        print "testSplitComplexes_thresholding"
        test_complexes, train_complexes, test_ppis, train_ppis, neg_test_ppis, neg_train_ppis = sc.split_complexes(self.complexes, size_threshold=2, threshold_fraction=0.0)
        print "test_complexes:"
        print test_complexes
        print "train_complexes:"
        print train_complexes
        print "test_ppis: "
        print test_ppis
        print "train_ppis: "
        print train_ppis
        print "neg_test_ppis: "
        print neg_test_ppis
        print "neg_train_ppis: "
        print neg_train_ppis
        print ""

        assert set(['i','g']) not in [set(x) for x in test_ppis]
        assert set(['a','i']) in [set(x) for x in test_ppis]

        for test_ppi in test_ppis:
            assert test_ppi not in [set(x) for train_complex in train_complexes for x in it.combinations(train_complex,2) ] 
        for train_ppi in train_ppis:
            assert train_ppi not in [set(x) for test_complex in test_complexes for x in it.combinations(test_complex,2) ]  

        test_complexes, train_complexes, test_ppis, train_ppis, neg_test_ppis, neg_train_ppis = sc.split_complexes(self.complexes, size_threshold=2, threshold_fraction=0.0, remove_large_complexes=True)
        print "remove large complexes"
        print "test_complexes:"
        print test_complexes
        print "train_complexes:"
        print train_complexes
        print "test_ppis: "
        print test_ppis
        print "train_ppis: "
        print train_ppis
        print "neg_test_ppis: "
        print neg_test_ppis
        print "neg_train_ppis: "
        print neg_train_ppis
        print ""

        for c in test_complexes:
            assert len(c) <= 2
        for c in train_complexes:
            assert len(c) <= 2

if __name__ == "__main__":
        unittest.main()


