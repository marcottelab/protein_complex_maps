#!/usr/bin/python

# Unit tests for shared bait feature (hypergeometric test)

import unittest
import protein_complex_maps.features.shared_bait_feature as sbf
import numpy as np
import pandas as pd


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

        print result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')]
        print result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count
        assert((result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].pair_count == 3).all())
        np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values, 2.302585 )



if __name__ == "__main__":
        unittest.main()


