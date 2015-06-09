#!/usr/bin/python

# Unit tests for purification recipe

import unittest
import protein_complex_maps.purification_recipe.purification_recipe as pr
import numpy as np
import pandas


class PurificationRecipeTest(unittest.TestCase):

	def setUp(self,):

            self.sample_pickle_filenames = ["./test_data.p", "./test_data2.p"]

        def testCreateObject(self,):

            proteins = ['ENSG00000010438','ENSG00000013275']
            percent_threshold = 0.9
            self.pr_obj = pr.PurificationRecipe( self.sample_pickle_filenames, proteins, percent_threshold )

            print "in testCreateObject"

            #kdrew: this is test_data.p data frame
            #ENSG00000008018  ENSG00000010438  ENSG00000013275
            #C101              500               50               10
            #C102              200              150                0
            #C103              300                0               90
            
            #kdrew: testing filename access
            assert(self.pr_obj.df_dict.keys()[1] == "./test_data.p")
            assert(self.pr_obj.df_dict.keys()[0] == "./test_data2.p")

            #kdrew: testing normalized data values, normalized so fractions sum to 1.0
            value = self.pr_obj.df_dict['./test_data.p']['ENSG00000010438'].values[0]
            value2 = self.pr_obj.df_dict['./test_data.p']['ENSG00000010438'].values[1]
            value3 = self.pr_obj.df_dict['./test_data.p']['ENSG00000010438'].values[2]
            np.testing.assert_almost_equal( value, 0.08928571428571429 )
            np.testing.assert_almost_equal( value2, 0.42857142857142855 )
            np.testing.assert_almost_equal( value3, 0.0 )

            passed_fractions = self.pr_obj.df_index_dict['./test_data.p']
            print passed_fractions
            #kdrew: only C101 has entries for both proteins, other fractions have 0.0 values for one or the other or both
            assert(passed_fractions[0] == 'C101')
            assert(len(passed_fractions) == 1)

        def testSingleProtein(self,):

            print "in testSingleProtein"
                
            proteins = ['ENSG00000010438']
            percent_threshold = 0.9
            self.pr_obj = pr.PurificationRecipe( self.sample_pickle_filenames, proteins, percent_threshold )
            self.pr_obj.create_recipes()

            #kdrew: testing calculated values
            np.testing.assert_almost_equal(self.pr_obj.results_list[4].purity_percent, 0.42857142857142855)
            np.testing.assert_almost_equal(self.pr_obj.results_list[4].yield_percent, 0.8275862068965518)
            assert(self.pr_obj.results_list[4].protein_percent, 1.0)
            assert(self.pr_obj.results_list[4].fractions[0] == 'C102')



            #self.pr_obj.show_results()

        def testMultiProtein(self,):

            print "in testDoubleProtein"

            proteins = ['ENSG00000010438','ENSG00000013275']
            percent_threshold = 0.9
            self.pr_obj = pr.PurificationRecipe( self.sample_pickle_filenames, proteins, percent_threshold )
            self.pr_obj.create_recipes()

            self.pr_obj.show_results()
                





if __name__ == "__main__":
	unittest.main()


