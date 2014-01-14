#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs
import numpy as np


class ReadDataTest(unittest.TestCase):

	def setUp(self,):

		sample_filename = "./test_data_uniprot.txt"
		sample_file = open(sample_filename, 'rb')
		self.data_matrix, self.name_list = rd.read_datafile(sample_file)

		print self.data_matrix

	def testRatio(self,):

		ratios = rs.calculate_ratio(self.data_matrix, self.name_list, "Q9Y6G5", "Q9UBI1")
		
		test_ratios = np.array([-1.01531676, -1.01531676, -0.8392255,  -1.14025549]) 
		print "test_ratios", test_ratios
		print ratios[0], test_ratios[0]
		#np.testing.assert_almost_equal( ratios[0], test_ratios[0] )  
		#np.testing.assert_almost_equal( ratios[1], test_ratios[1] )  
		#np.testing.assert_almost_equal( ratios[2], test_ratios[2] )  
		#np.testing.assert_almost_equal( ratios[3], test_ratios[3] )  
		np.testing.assert_almost_equal( ratios, test_ratios )  
		

		#assert(c_dmat[np.ix_([index2],[7])] == 0.0)


if __name__ == "__main__":
	unittest.main()


