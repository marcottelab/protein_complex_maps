#!/usr/bin/python

# Unit tests for purification recipe

import unittest
import protein_complex_maps.purification_recipe.purification_recipe as pr
import numpy as np
import pandas


class PurificationRecipeTest(unittest.TestCase):

	def setUp(self,):

            self.sample_pickle_filenames = ["./test_data.p", "./test_data2.p"]

            self.sample_reorder_pickle_filenames = ["./test_data.p", "./test_data2_reorder_proteins.p"]

            self.sample_uniprot_pickle_filenames = ["./test_data_ids_mapped.p", "./test_data2_ids_mapped.p"]

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
		assert(self.pr_obj.results_list[4].protein_percent == 1.0)
		assert(self.pr_obj.results_list[4].fractions[0] == 'C102')



		#self.pr_obj.show_results()

	def testMultiProtein(self,):

		print "in testDoubleProtein"

		proteins = ['ENSG00000010438','ENSG00000013275']
		percent_threshold = 0.9
		#kdrew: not normalizing so it is easier to manually calculate results
		self.pr_obj = pr.PurificationRecipe( self.sample_pickle_filenames, proteins, percent_threshold, normalize_by_fraction_flag=False )
		self.pr_obj.create_recipes()

		self.pr_obj.show_results()

		#print self.pr_obj.results_list[3]

		##ProtID TotalCount  C101    C102    C103    
		#ENSG00000008018 1000.0  500.0   200.0   300.0
		#ENSG00000010438 200.0   50.0    150.0   0.0 
		#ENSG00000013275 100.0   10.0    0.0 90.0    
		
		##ProtID TotalCount  C201    C202    C203    
		#ENSG00000008018 1000.0  100.0   0.0 900.0
		#ENSG00000010438 200.0   50.0    100.0   50.0    
		#ENSG00000013275 100.0   50.0   50.0 0.0 

		#('C201', 'C101')
		#ENSG00000008018    50.0
		#ENSG00000010438    12.5
		#ENSG00000013275     5.0


		#kdrew: C201 -> C101 
		#kdrew: for each protein multiply probability of seeing protein in fraction, 
		#kdrew: represents the expected number of proteins after a series of columns
		#kdrew: P(prot|C201) * P(prot|C101) * Total_Protein
		#kdrew: 100.0/1000 * 500.0/1000 = 0.05
		#kdrew: 1000 * 0.05 = 50.0
		#kdrew: 50.0/200 * 50.0/200 = 0.0625
		#kdrew: 200 * 0.0625 = 12.5
		#kdrew: 50.0/100 * 10.0/100 = 0.05 
		#kdrew: 100 * 0.05a = 5.0
		
		#kdrew: purity_percent:  (5.0 + 12.5) / (5.0 + 12.5 + 50.0) = 0.25925925925925924

		#kdrew: yield = mean( P(prot|C201) * P(prot|C101)) 
		#kdrew: yield 12.5/200, 5.0/100

		np.testing.assert_almost_equal(self.pr_obj.results_list[3].purity_percent, 0.25925925925925924)
		assert(self.pr_obj.results_list[3].protein_percent == 1.0)
		np.testing.assert_almost_equal(self.pr_obj.results_list[3].yield_percent, 0.056250000000000001)


	def testMultiProtein_reorder(self,):

		print "in testReorderProtein"

		proteins = ['ENSG00000010438','ENSG00000013275']
		percent_threshold = 0.9
		#kdrew: not normalizing so it is easier to manually calculate results
		self.pr_obj = pr.PurificationRecipe( self.sample_reorder_pickle_filenames, proteins, percent_threshold, normalize_by_fraction_flag=False )
		self.pr_obj.create_recipes()

		self.pr_obj.show_results()

		np.testing.assert_almost_equal(self.pr_obj.results_list[3].purity_percent, 0.25925925925925924)
		assert(self.pr_obj.results_list[3].protein_percent == 1.0)
		np.testing.assert_almost_equal(self.pr_obj.results_list[3].yield_percent, 0.056250000000000001)


	def testMultiProtein_uniprot(self,):

		print "in testUniprot"

		#proteins = ['ENSG00000010438','ENSG00000013275']
		proteins = ['P35030','P43686']
		percent_threshold = 0.9
		#kdrew: not normalizing so it is easier to manually calculate results
		self.pr_obj = pr.PurificationRecipe( self.sample_uniprot_pickle_filenames, proteins, percent_threshold, normalize_by_fraction_flag=False )
		self.pr_obj.create_recipes()

		self.pr_obj.show_results()

		np.testing.assert_almost_equal(self.pr_obj.results_list[3].purity_percent, 0.25925925925925924)
		assert(self.pr_obj.results_list[3].protein_percent == 1.0)
		np.testing.assert_almost_equal(self.pr_obj.results_list[3].yield_percent, 0.056250000000000001)

	def testMultiProtein_uniprot_normalized(self,):

		print "in testUniprot_normalized"

		#proteins = ['ENSG00000010438','ENSG00000013275']
		proteins = ['P35030','P43686']
		percent_threshold = 0.9
		#kdrew: testing normalized 
		self.pr_obj = pr.PurificationRecipe( self.sample_uniprot_pickle_filenames, proteins, percent_threshold, normalize_by_fraction_flag=True )
		self.pr_obj.create_recipes()

		self.pr_obj.show_results()

		#kdrew: not actually testing anything, need to compute results by hand
		#np.testing.assert_almost_equal(self.pr_obj.results_list[3].purity_percent, 0.25925925925925924)
		#assert(self.pr_obj.results_list[3].protein_percent == 1.0)
		#np.testing.assert_almost_equal(self.pr_obj.results_list[3].yield_percent, 0.056250000000000001)


if __name__ == "__main__":
	unittest.main()


