#!/usr/bin/python

# Unit tests for stoichiometry functions

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.protein_util as pu
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs
import protein_complex_maps.stoichiometry.stoichiometry as st
import numpy as np


class StoichiometryTest(unittest.TestCase):

	def setUp(self,):

		filename = "./stoichiometry_list.txt"
		sample_file = open(filename, 'rb')
		self.s_set = st.Stoichiometries()
		self.s_set.read_stoichiometries(sample_file)

		#print self.data_matrix

	def testStoichiometrySet(self,):

		for i in self.s_set:
			print i, i.count

		assert(len(self.s_set) == 20)


if __name__ == "__main__":
	unittest.main()


