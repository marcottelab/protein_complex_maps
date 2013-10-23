#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.bicluster.bicluster as bc
import numpy as np


class BiclusterTest(unittest.TestCase):

	def setUp(self,):

		self.data_matrix = np.arange(50).reshape(10, 5)
		#print self.data_matrix
		self.bicluster1 = bc.Bicluster()
		self.bicluster2 = bc.Bicluster(rows=[0,1,7,8],cols=[0,3])

		self.bicluster1.add_row(1)
		self.bicluster1.add_column(2)

	def testAdd(self, ):
		assert( self.bicluster1.rows() == [1] )
		assert( self.bicluster1.columns() == [2] )


	def testSubMatrix(self,):
		submat = self.bicluster1.get_submatrix(self.data_matrix)
		assert( submat == np.array([7]).reshape(1,1))

		submat = self.bicluster2.get_submatrix(self.data_matrix)
		#print submat
		submat_comp = np.array([0,3,5,8,35,38,40,43]).reshape(4,2)
		#print submat_comp
		assert( submat[0,0] == submat_comp[0,0] )
		assert( submat[0,1] == submat_comp[0,1] )
		assert( submat[1,0] == submat_comp[1,0] )
		assert( submat[1,1] == submat_comp[1,1] )
		assert( submat[2,0] == submat_comp[2,0] )
		assert( submat[2,1] == submat_comp[2,1] )
		assert( submat[3,0] == submat_comp[3,0] )
		assert( submat[3,1] == submat_comp[3,1] )

	def testSubMatrixWO(self,):
		submat = self.bicluster2.get_submatrix(self.data_matrix, without_rows=[1])
		#print submat
		submat_comp = np.array([0,3,5,8,35,38,40,43]).reshape(4,2)
		#print submat_comp
		assert( submat[0,0] == submat_comp[0,0] )
		assert( submat[0,1] == submat_comp[0,1] )
		assert( submat[1,0] == submat_comp[2,0] )
		assert( submat[1,1] == submat_comp[2,1] )
		assert( submat[2,0] == submat_comp[3,0] )
		assert( submat[2,1] == submat_comp[3,1] )

		submat = self.bicluster2.get_submatrix(self.data_matrix, without_cols=[0])
		#print submat
		submat_comp = np.array([0,3,5,8,35,38,40,43]).reshape(4,2)
		#print submat_comp
		assert( submat[0,0] == submat_comp[0,1] )
		assert( submat[1,0] == submat_comp[1,1] )
		assert( submat[2,0] == submat_comp[2,1] )
		assert( submat[3,0] == submat_comp[3,1] )

	def testOutsideRowsCols(self,):
		out_rows = self.bicluster2.get_outside_rows(self.data_matrix)
		#print "out_rows", out_rows
		assert(out_rows == set([2,3,4,5,6,9]))

		out_cols = self.bicluster2.get_outside_cols(self.data_matrix)
		assert(out_cols == set([1,2,4]))

	def testRandOutsideRowsCols(self,):
		rand_rows = self.bicluster2.get_random_outside_rows(self.data_matrix,seed=0)
		assert( rand_rows == [2,3,5,9] )

		rand_cols = self.bicluster2.get_random_outside_cols(self.data_matrix,seed=0)
		assert( rand_cols == [2,4] )

	def testScale(self,):
		print self.bicluster2
		print self.data_matrix
		dmat = self.bicluster2.scale(self.data_matrix,100)
		print dmat

		assert( dmat[0,0] == 0 )
		assert( dmat[1,0] == 500 )
		assert( dmat[2,0] == 10 )
		assert( dmat[7,0] == 3500 )
		assert( dmat[8,0] == 4000 )
		assert( dmat[0,3] == 300 )
		assert( dmat[1,3] == 800 )
		assert( dmat[2,3] == 13 )
		assert( dmat[7,3] == 3800 )
		assert( dmat[8,3] == 4300 )
		
		

if __name__ == "__main__":
	unittest.main()


