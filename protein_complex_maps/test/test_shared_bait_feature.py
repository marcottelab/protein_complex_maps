#!/usr/bin/python

# Unit tests for comparing clusters to gold standard complexes

import unittest
import protein_complex_maps.features.shared_bait_feature as sbf
import numpy as np
import pandas as pd
import scipy


class SharedBaitFeatureTest(unittest.TestCase):

    def setUp(self,):

        self.feature_dict = dict()
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
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'B') & (result_table['gene_id2'] == '2')].neg_ln_pval.values[0], 0.510826, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '5') & (result_table['gene_id2'] == 'C')].neg_ln_pval.values[0], 1.609438, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'E') & (result_table['gene_id2'] == '8')].pval.values[0], 0.4 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == 'A')].pair_count.values[0], 1 )

        #print result_table.query("gene_id1 == '2' and gene_id2 == '3'")
        np.testing.assert_almost_equal( result_table.query("gene_id1 == '2' and gene_id2 == '3'").neg_ln_pval.values, 2.302585 )

    def testSharedBaitFeatureExponents(self, ):
        self.feature_dict['gene_id'] = [1,2,2,3,3,4,5,6,2,3,7,8,8,'A','B','C','D','E','F','2','G','3']
        self.feature_dict['bait_geneid'] = ['B','B','A','A','B','C','C','C','C','C','D','D','E','A','B','C','D','E','F','F','G','G']
        self.feature_table = pd.DataFrame(self.feature_dict)
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id', logchoose=True)
        print "result_table logchoose"
        print result_table
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id', transform_pvalue=True)
        print "result_table exponents"
        print result_table
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].exponents.values[0][0], -1.0704414117, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].exponents.values[0][1], -3.55534806149, decimal=6 )

        min_exponent = min(result_table.exponents.apply(min))
        print "min_exponent: %s" % min_exponent
        result_table = result_table.exponents.apply(lambda x: [xx-min_exponent for xx in x])
        print result_table

    def testSharedBaitFeatureLogChoose(self, ):
        result_table = sbf.shared_bait_feature(self.feature_table, 'bait_geneid', 'gene_id', logchoose=True)
        print "result_table exponents"
        print result_table
        assert((result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count == 3).all())
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values[0], 2.302585, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'B') & (result_table['gene_id2'] == '2')].neg_ln_pval.values[0], 0.510826, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '5') & (result_table['gene_id2'] == 'C')].neg_ln_pval.values[0], 1.609438, decimal=6 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == 'E') & (result_table['gene_id2'] == '8')].pval.values[0], 0.4 )
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == 'A')].pair_count.values[0], 1 )

    def testSharedBaitFeatureBoundaryCase(self, ):

        print "Boundary Case"
        print "k=670,m=1196,n=813,N=5891"
        pval = sbf.pval(k=670,m=1196,n=813,N=5891, logchoose=True)
        print pval
        exponents = sbf.hypergeometric_exponents(k=670,m=1196,n=813,N=5891)
        print exponents
        max_exp = np.max(exponents)
        print max_exp
        transformed_exponents = [ex-max_exp for ex in exponents]
        print "transformed_exponents"
        print transformed_exponents
        transformed_pval = sum([scipy.exp(ex) for ex in transformed_exponents])
        print "transformed pval"
        print transformed_pval
        pval2 = transformed_pval * scipy.exp(max_exp)
        print "pval2"
        print pval2
        print sbf.transform_pval(670,1196,813,5891,670)
        print "transformed pval"
        print sbf.transform_pval(670,1196,813,5891,-500, calc_true_pval=False)

        

        print "k=722,m=1088,n=1064,N=5891"
        pval = sbf.pval(k=722,m=1088,n=1064,N=5891, logchoose=True)
        print pval
        exponents= sbf.hypergeometric_exponents(k=722,m=1088,n=1064,N=5891)
        print exponents
        print np.max(exponents)

        print "k=22,m=1088,n=1064,N=5891"
        pval = sbf.pval(k=22,m=1088,n=1064,N=5891, logchoose=True)
        print pval
        exponents= sbf.hypergeometric_exponents(k=22,m=1088,n=1064,N=5891)
        print exponents
        print np.max(exponents)
        print "transformed pval"
        print sbf.transform_pval(22,1088,1064,5891,-100, calc_true_pval=False)

        print "k=3,m=5,n=4,N=5891"
        pval = sbf.pval(k=3,m=5,n=4,N=5891, logchoose=True)
        print pval
        exponents= sbf.hypergeometric_exponents(k=3,m=5,n=4,N=5891)
        print exponents
        print np.max(exponents)
        print "transformed pval"
        print sbf.transform_pval(3,5,4,5891,-100, calc_true_pval=False)

    def testSharedBaitFeatureTransformPval(self, ):
        pval = sbf.pval(k=10,m=12,n=11,N=100, logchoose=False)
        logpval = sbf.pval(k=10,m=12,n=11,N=100, logchoose=True)
        transpval = sbf.transform_pval(k=10,m=12,n=11,N=100,transform_val=10)

        print "pval: %s, logpval: %s, transpval: %s" % (pval, logpval, transpval)
        
        np.testing.assert_almost_equal(pval, logpval)
        np.testing.assert_almost_equal(pval, transpval)
        np.testing.assert_almost_equal(logpval, transpval)

if __name__ == "__main__":
        unittest.main()


