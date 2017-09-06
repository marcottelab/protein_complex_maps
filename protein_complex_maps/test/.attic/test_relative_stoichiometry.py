#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.protein_util as pu
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs
import protein_complex_maps.stoichiometry.stoichiometry as st
import numpy as np


class ReadDataTest(unittest.TestCase):

	def setUp(self,):

		sample_filename = "./test_data_uniprot.txt"
		sample_file = open(sample_filename, 'rb')
		self.msds = rd.MSDataSet()
		self.msds.load_file(sample_file, header=True)

		#print self.data_matrix

	def testRatio(self,):

		#p1_length = pu.get_length_uniprot("Q9Y6G5")
		#p2_length = pu.get_length_uniprot("Q9UBI1")

		dm = nu.normalize_length( self.msds.get_data_matrix(), self.msds.get_id_dict() )
		self.msds.set_data_matrix(dm)
		ratios = rs.calculate_ratio(self.msds, "Q9Y6G5", "Q9UBI1" )
		
		test_ratios = np.array([-1.01531676, -1.01531676, -0.8392255,  -1.14025549]) 
		print "test_ratios", test_ratios
		print ratios[0], test_ratios[0]
		#np.testing.assert_almost_equal( ratios[0], test_ratios[0] )  
		#np.testing.assert_almost_equal( ratios[1], test_ratios[1] )  
		#np.testing.assert_almost_equal( ratios[2], test_ratios[2] )  
		#np.testing.assert_almost_equal( ratios[3], test_ratios[3] )  
		np.testing.assert_almost_equal( ratios, test_ratios )  
		
		#assert(c_dmat[np.ix_([index2],[7])] == 0.0)
	
	def testRelativeStoichiometryProbability(self,):

		stoichiometry = dict()
		stoichiometry['A'] = 2 
		stoichiometry['B'] = 20

		prior = 0.2

		prot_ids = dict()
		prot_ids['A'] = "Q9Y6G5"
		prot_ids['B'] = "Q9UBI1"

		#kdrew: un-normalized length
		lp = rs.relative_stoichiometry_probability(stoichiometry, prior, self.msds, prot_ids, mean_ratio=False, median_ratio=False, set_std=False)
		np.testing.assert_almost_equal( lp[0], -5.30850095496 )
		np.testing.assert_almost_equal( lp[1], 4 )

	def testRelativeStoichiometryProbability2(self,):

		stoichiometry = dict()
		stoichiometry['A'] = 2 
		stoichiometry['B'] = 20

		prior = 0.2

		prot_ids = dict()
		prot_ids['A'] = "Q9Y6G5"
		prot_ids['B'] = "Q9UBI1"

		#kdrew: normalized length
		dm = nu.normalize_length( self.msds.get_data_matrix(), self.msds.get_id_dict() )
		self.msds.set_data_matrix(dm)

		lp = rs.relative_stoichiometry_probability(stoichiometry, prior, self.msds, prot_ids, mean_ratio=False,median_ratio=False,set_std=False)
		np.testing.assert_almost_equal( lp[0], -5.30818667031 )
		np.testing.assert_almost_equal( lp[1], 4 )


	def testRelativeStoichiometry(self,):
		print "relative_stoichiometry"
		prot_ids = ["Q9Y6G5", "Q9UBI1"]

		filename = "./stoichiometry_list.txt"
		sample_file = open(filename, 'rb')
		s_set = st.Stoichiometries() 
		s_set.read_stoichiometries(sample_file)

		results = rs.relative_stoichiometry(self.msds, prot_ids, s_set)
		print results


if __name__ == "__main__":
	unittest.main()


