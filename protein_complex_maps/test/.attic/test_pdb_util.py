#!/usr/bin/python

# Unit tests for pdb util 

import unittest
import protein_complex_maps.pdb_util as pu
import numpy as np
import Bio.PDB


class PDBUtilTest(unittest.TestCase):

	def setUp(self,):

		#kdrew: parse filename to get pdbid
		pdb_file = "./1b70.pdb"
		pdbid = "1b70"
		self.chain1 = 'A'
		self.chain2 = 'B'
		try:
			self.structure = Bio.PDB.PDBParser().get_structure(pdbid, pdb_file)
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
			return

		try:
			self.chain_one = self.structure[0][self.chain1]
		except KeyError:
			print "missing chain: %s" % (self.chain1,)
			return

		try:
			self.chain_two = self.structure[0][self.chain2]
		except KeyError:
			print "missing chain: %s" % (self.chain2,)
			return


	def testResidueDistance(self,):
		r1 = self.chain_one[86]
		r2 = self.chain_two[1]
		chain_Nterm_dist = pu.calc_residue_dist(r1, r2)
		print chain_Nterm_dist
		np.testing.assert_almost_equal( chain_Nterm_dist, 76.49939, decimal=5 )

	def testChainDistance(self,):
		chain_dist_mat = pu.calc_dist_matrix(self.chain_one, self.chain_two)
		print chain_dist_mat[0][0]
		np.testing.assert_almost_equal( chain_dist_mat[0][0], 76.49939, decimal=5 )

	def testChainDistanceComplete(self,):
		chain_dist_mat = pu.calc_dist_matrix_complete(self.chain_one, self.chain_two)
		print chain_dist_mat[0][265]
		np.testing.assert_almost_equal( chain_dist_mat[0][265], 76.49939, decimal=5 )

if __name__ == "__main__":
	unittest.main()


