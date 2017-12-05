#!/usr/bin/python

# Unit tests for shift_frac

import unittest
import protein_complex_maps.features.ExtractFeatures.Features as eff
import protein_complex_maps.features.shift_frac as sf
import numpy as np
import pandas as pd


class ShiftFracTest(unittest.TestCase):

    def setUp(self,):

        d = dict()
        d['ids'] = ['a','b','c','d','e']
        d['frac1'] = [1.0,2.0,0.0,0.0,0.0]
        d['frac2'] = [0.0,0.0,3.0,4.0,0.0]
        d['frac3'] = [2.0,1.0,0.0,0.0,7.0]
        d['frac4'] = [0.0,0.0,6.0,8.0,10.0]
        df = pd.DataFrame(d, dtype=float)
        df = df.set_index('ids')
        self.elution = eff.Elut()
        self.elution.df = df

        d2 = dict()
        d2['ids'] = ['a','b','c','d']
        d2['frac1'] = [1.0,1.0,0.0,2.0]
        d2['frac2'] = [0.0,0.0,1.0,4.0]
        d2['frac3'] = [2.0,1.0,0.0,0.0]
        d2['frac4'] = [0.0,0.0,6.0,10.0]
        df2 = pd.DataFrame(d2, dtype=float)
        df2 = df2.set_index('ids')
        self.elution2 = eff.Elut()
        self.elution2.df = df2

    def testShiftFrac(self, ):
        #kdrew: subtracts one dataframe from another and sums the absolute values
        shift_frac_sum = sf.calc_shift_frac(self.elution, self.elution2)
        print shift_frac_sum
        assert(shift_frac_sum.loc['a'] == 0.0)
        assert(shift_frac_sum.loc['b'] == 1.0)
        assert(shift_frac_sum.loc['c'] == 2.0)
        assert(shift_frac_sum.loc['d'] == 4.0)
        assert(shift_frac_sum.loc['e'] == 17.0)

        #np.testing.assert_almost_equal(result_table[(result_table['gene_id1'] == '2') & (result_table['gene_id2'] == '3')].neg_ln_pval.values, 2.302585 )



if __name__ == "__main__":
        unittest.main()


