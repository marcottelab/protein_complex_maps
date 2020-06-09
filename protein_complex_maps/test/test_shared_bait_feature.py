#!/usr/bin/python

# Unit tests for shared bait feature (hypergeometric test)

import unittest
import protein_complex_maps.features.shared_bait_feature as sbf
import numpy as np
import pandas as pd
import scipy


class SharedBaitFeatureTest(unittest.TestCase):

    def setUp(self,):

        self.feature_dict = dict()
        #self.feature_dict['abundance'] = [10.1, 90.2, 2.5, 50.7, 88.8, 45.0, 22.1, 44.3, 76.9, 93.2, 67.3, 21.1, 14.4, 59.0, 88.7, 78.7, 85.5, 95.7]
        self.feature_dict['abundance'] = [1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 1.0, 1.0, 6.0, 6.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        self.feature_dict['abundance_present'] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        self.feature_dict['gene_id'] = [1,2,2,3,3,4,5,6,2,3,7,8,8,'A','B','C','D','E']
        self.feature_dict['bait_geneid'] = ['B','B','A','A','B','C','C','C','C','C','D','D','E','A','B','C','D','E']
        self.feature_table = pd.DataFrame(self.feature_dict)

    def testSharedBaitFeature(self, ):
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id')
        print "result_table"
        print result_table

        #print "result_table 2,3"
        #print result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')]
        #print "result_table 2,3 pair_count"
        #print result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count
        assert((result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count == 3).all())
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values[0], 2.302585, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == 'B')].neg_ln_pval.values[0], 0.510826, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '5') & (result_table['gene_id2'] == 'C')].neg_ln_pval.values[0], 1.609438, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '8') & (result_table['gene_id2'] == 'E')].pval.values[0], 0.4 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'A') & (result_table['gene_id2'] == '2')].pair_count.values[0], 1 )

        #print result_table.query("gene_id1 == '2' and gene_id2 == '3'")
        np.testing.assert_almost_equal( result_table.query("gene_id1 == '2' and gene_id2 == '3'").neg_ln_pval.values, 2.302585 )

    def testSharedBaitFeatureUseAbundancePresent(self, ):
        #kdrew: tests the abundance code using a single count per experiment, this should result in the same p-values as using presence absence 
        print "testSharedBaitFeatureUseAbundancePresent"
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id', 'abundance_present',use_abundance=True)
        print "result_table"
        print result_table

        assert((result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count == 3).all())
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values[0], 2.302585, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == 'B')].neg_ln_pval.values[0], 0.510826, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '5') & (result_table['gene_id2'] == 'C')].neg_ln_pval.values[0], 1.609438, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '8') & (result_table['gene_id2'] == 'E')].pval.values[0], 0.4 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'A') & (result_table['gene_id2'] == '2')].pair_count.values[0], 1 )

        #print result_table.query("gene_id1 == '2' and gene_id2 == '3'")
        np.testing.assert_almost_equal( result_table.query("gene_id1 == '2' and gene_id2 == '3'").neg_ln_pval.values, 2.302585 )
 
    def testSharedBaitFeatureUseAbundance(self, ):
        #kdrew: tests abundance code with pseudo real abundances, 
        #kdrew: should fail test currently
        print "testSharedBaitFeatureUseAbundance"
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id', 'abundance',use_abundance=True)
        print "result_table"
        print result_table

        assert((result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count == 11).all())
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pval.values[0], 0.11029411764705878, decimal=6 )
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'B') & (result_table['gene_id2'] == '2')].neg_ln_pval.values[0], 0.510826, decimal=6 )
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '5') & (result_table['gene_id2'] == 'C')].neg_ln_pval.values[0], 1.609438, decimal=6 )
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'E') & (result_table['gene_id2'] == '8')].pval.values[0], 0.4 )
        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == 'A')].pair_count.values[0], 1 )

        #print result_table.query("gene_id1 == '2' and gene_id2 == '3'")
        np.testing.assert_almost_equal( result_table.query("gene_id1 == '2' and gene_id2 == '3'").pval.values, 0.11029411764705878)

    def testSharedBaitFeatureBoundaryCaseOld(self, ):

        print "OLD Boundary Case"
        print "k=670,m=1196,n=813,N=5891"
        pval = sbf.pval_old(k=670,m=1196,n=813,N=5891, logchoose=True)
        print "*****************"
        print pval

        print "k=722,m=1088,n=1064,N=5891"
        pval = sbf.pval_old(k=722,m=1088,n=1064,N=5891, logchoose=True)
        print pval

        print "k=22,m=1088,n=1064,N=5891"
        pval = sbf.pval_old(k=22,m=1088,n=1064,N=5891, logchoose=True)
        print pval

        print "k=3,m=5,n=4,N=5891"
        pval = sbf.pval_old(k=3,m=5,n=4,N=5891, logchoose=True)
        print pval

    def testSharedBaitFeatureBoundaryCase(self, ):

        print "Boundary Case"
        print "k=670,m=1196,n=813,N=5891"
        pval = sbf.pval(k=670,m=1196,n=813,N=5891)
        print "*****************"
        print pval
        #kdrew: integration test rather than unit test (hard to verify correct value without arbitrary machine precision
        np.testing.assert_almost_equal(pval, 1.56133062575003e-394)

        print "k=722,m=1088,n=1064,N=5891"
        pval = sbf.pval(k=722,m=1088,n=1064,N=5891)
        print pval

        print "k=22,m=1088,n=1064,N=5891"
        pval = sbf.pval(k=22,m=1088,n=1064,N=5891)
        print pval

        print "k=3,m=5,n=4,N=5891"
        pval = sbf.pval(k=3,m=5,n=4,N=5891)
        print pval
        np.testing.assert_almost_equal(pval, 1.17423417751186e-9)

    def testSharedBaitFeatureTransformPvalOld(self, ):
        pval = sbf.pval_old(k=10,m=12,n=11,N=100, logchoose=False)
        
        np.testing.assert_almost_equal(pval, 4.1093045454984257e-11)

    def testSharedBaitFeaturePval(self, ):
        pval = sbf.pval(k=10,m=12,n=11,N=100)
        np.testing.assert_almost_equal(pval, 4.1093045454984257e-11)

    def testSharedBaitFeaturePval_abundance(self, ):
        pval = sbf.pval(k=10,m=12,n=11,N=100)
        np.testing.assert_almost_equal(pval, 4.1093045454984257e-11)


if __name__ == "__main__":
        unittest.main()


