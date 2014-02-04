#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import numpy as np


class NormalizationTest(unittest.TestCase):

	def setUp(self,):

		self.data_matrix = np.matrix( np.zeros((10,10)) )

		self.data_matrix[4,5] = 1.9
		self.data_matrix[4,8] = 2.9
		self.data_matrix[2,5] = 3.9
		self.data_matrix[2,8] = 4.9
		self.data_matrix[2,6] = 5.9
		self.data_matrix[7,8] = 6.9

		#print self.data_matrix
		#kdrew: original matrix
		#[[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   3.9  5.9  0.   4.9  0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   1.9  0.   0.   2.9  0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   6.9  0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]]


	def testZero(self, ):
		data_matrix_zero_row = nu.remove_zero(self.data_matrix, zero_rows=True, zero_columns=False)
		#print data_matrix_zero_row
		#kdrew: removed rows with all zeros
		#[[ 0.   0.   0.   0.   0.   3.9  5.9  0.   4.9  0. ]
	    #[ 0.   0.   0.   0.   0.   1.9  0.   0.   2.9  0. ]
		#[ 0.   0.   0.   0.   0.   0.   0.   0.   6.9  0. ]]
		assert( data_matrix_zero_row[0,5] == 3.9 )
		assert( data_matrix_zero_row[0,6] == 5.9 )
		assert( data_matrix_zero_row[0,8] == 4.9 )
		assert( data_matrix_zero_row[2,8] == 6.9 )
		assert( data_matrix_zero_row[1,5] == 1.9 )
		assert( data_matrix_zero_row[1,8] == 2.9 )

		data_matrix_zero_column = nu.remove_zero(self.data_matrix, zero_rows=False, zero_columns=True)
		#print data_matrix_zero_column
		#kdrew: removed columns with all zeros
		#[[ 0.   0.   0. ]
		#[ 0.   0.   0. ]
		#[ 3.9  5.9  4.9]
		#[ 0.   0.   0. ]
		#[ 1.9  0.   2.9]
		#[ 0.   0.   0. ]
		#[ 0.   0.   0. ]
		#[ 0.   0.   6.9]
		#[ 0.   0.   0. ]
		#[ 0.   0.   0. ]]
		assert( data_matrix_zero_column[2,0] == 3.9 )
		assert( data_matrix_zero_column[2,1] == 5.9 )
		assert( data_matrix_zero_column[2,2] == 4.9 )
		assert( data_matrix_zero_column[4,0] == 1.9 )
		assert( data_matrix_zero_column[4,2] == 2.9 )
		assert( data_matrix_zero_column[7,2] == 6.9 )

		data_matrix_zero = nu.remove_zero(self.data_matrix, zero_rows=True, zero_columns=True)
		#print data_matrix_zero
		#kdrew: removed columns and rows that had all zeros
		#[[ 3.9  5.9  4.9]
		#[ 1.9  0.   2.9]
		#[ 0.   0.   6.9]]
		assert( data_matrix_zero[0,0] == 3.9 )
		assert( data_matrix_zero[0,1] == 5.9 )
		assert( data_matrix_zero[0,2] == 4.9 )
		assert( data_matrix_zero[1,0] == 1.9 )
		assert( data_matrix_zero[1,2] == 2.9 )
		assert( data_matrix_zero[2,2] == 6.9 )


	def testNormalize(self,):
		#print self.data_matrix
		data_matrix_normalized_rows = nu.normalize_over_rows(self.data_matrix)
		#print data_matrix_normalized_rows
		#kdrew: normalized where individual columns sum to 1.0
		#[[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.67241379 1.          0.          0.33333333  0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.32758621 0.          0.          0.19727891  0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.46938776  0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]]
		np.testing.assert_almost_equal( data_matrix_normalized_rows[2,5], 0.67241379 )
		np.testing.assert_almost_equal( data_matrix_normalized_rows[2,6], 1.0 )
		np.testing.assert_almost_equal( data_matrix_normalized_rows[2,8], 0.33333333 )
		np.testing.assert_almost_equal( data_matrix_normalized_rows[4,5], 0.32758621 )
		np.testing.assert_almost_equal( data_matrix_normalized_rows[4,8], 0.19727891 )
		np.testing.assert_almost_equal( data_matrix_normalized_rows[7,8], 0.46938776 )

		#print self.data_matrix
		data_matrix_normalized_columns = nu.normalize_over_columns(self.data_matrix)
		#print data_matrix_normalized_columns
		#kdrew: normalized where individual rows sum to 1.0
		#[[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.26530612 0.40136054  0.          0.33333333  0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.39583333 0.          0.          0.60416667  0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          1.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]
		#[ 0.          0.          0.          0.          0.          0.          0.  0.          0.          0.        ]]
		np.testing.assert_almost_equal( data_matrix_normalized_columns[2,5], 0.26530612 )
		np.testing.assert_almost_equal( data_matrix_normalized_columns[2,6], 0.40136054 )
		np.testing.assert_almost_equal( data_matrix_normalized_columns[2,8], 0.33333333 )
		np.testing.assert_almost_equal( data_matrix_normalized_columns[4,5], 0.39583333 )
		np.testing.assert_almost_equal( data_matrix_normalized_columns[4,8], 0.60416667 )
		np.testing.assert_almost_equal( data_matrix_normalized_columns[7,8], 1.0 )

	def testNoise(self,):
		data_matrix_column_noise = nu.add_noise_over_columns(self.data_matrix)
		#print data_matrix_column_noise
		#kdrew: added noise 1/(number of columns)
		#[[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  4.   6.   0.1  5.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  2.   0.1  0.1  3.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  7.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]]
		assert( data_matrix_column_noise[0,0] == 0.1 )
		assert( data_matrix_column_noise[2,5] == 4.0 )
		assert( data_matrix_column_noise[2,6] == 6.0 )
		assert( data_matrix_column_noise[2,7] == 0.1 )
		assert( data_matrix_column_noise[2,8] == 5.0 )
		assert( data_matrix_column_noise[4,5] == 2.0 )
		assert( data_matrix_column_noise[4,8] == 3.0 )
		assert( data_matrix_column_noise[7,8] == 7.0 )

		data_matrix_column_noise = nu.add_noise_over_rows(self.data_matrix)
		#print data_matrix_column_noise
		#kdrew: added noise 1/(number of rows), looks the same as above because both number of rows and columns is 10
		#[[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  4.   6.   0.1  5.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  2.   0.1  0.1  3.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  7.   0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]
		#[ 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1]]
		assert( data_matrix_column_noise[0,0] == 0.1 )
		assert( data_matrix_column_noise[2,5] == 4.0 )
		assert( data_matrix_column_noise[2,6] == 6.0 )
		assert( data_matrix_column_noise[2,7] == 0.1 )
		assert( data_matrix_column_noise[2,8] == 5.0 )
		assert( data_matrix_column_noise[4,5] == 2.0 )
		assert( data_matrix_column_noise[4,8] == 3.0 )
		assert( data_matrix_column_noise[7,8] == 7.0 )


	def testMinToOne(self,):
		data_matrix_column_noise = nu.add_noise_over_columns(self.data_matrix)
		print data_matrix_column_noise
		data_matrix_normalized = nu.normalize_over_rows(data_matrix_column_noise)
		print data_matrix_normalized
		data_matrix_mintoone = nu.min_to_one_scale(data_matrix_normalized)
		print data_matrix_mintoone
		print data_matrix_mintoone[0,0]
		np.testing.assert_almost_equal( data_matrix_mintoone[0,0], 15.7 )
		np.testing.assert_almost_equal( data_matrix_mintoone[1,8], 1.0 )
		np.testing.assert_almost_equal( data_matrix_mintoone[2,8], 50.0 )

	def testPoisson(self,):
		np.random.seed(124)
		data_matrix_poisson = nu.sample_noise(self.data_matrix, sample_module=np.random.poisson)		
		#print data_matrix_poisson
		#kdrew: added poisson noise
		#[[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 4 6 0 4 0]
		#[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 4 0 0 2 0]
		#[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 0 0 0 9 0]
		#[0 0 0 0 0 0 0 0 0 0]
		#[0 0 0 0 0 0 0 0 0 0]]
		assert( data_matrix_poisson[0,0] == 0.0 )
		assert( data_matrix_poisson[2,5] == 4.0 )
		assert( data_matrix_poisson[2,6] == 6.0 )
		assert( data_matrix_poisson[2,7] == 0.0 )
		assert( data_matrix_poisson[2,8] == 4.0 )
		assert( data_matrix_poisson[4,5] == 4.0 )
		assert( data_matrix_poisson[4,8] == 2.0 )
		assert( data_matrix_poisson[7,8] == 9.0 )

	def testLength(self,):
		
		id_dict = dict()
		id_dict['Q9Y6G5'] = 2
		id_dict['Q9UBI1'] = 4
		dm_length = nu.normalize_length(self.data_matrix, id_dict)
		print dm_length
		print dm_length[0,0]
		assert( np.isnan(dm_length[0,0]) )
		np.testing.assert_almost_equal( dm_length[2,5], 0.01930693 )

	def testPeptide(self, ):

		sample_filename2 = "./test_data2.txt"
		sample_file2 = open(sample_filename2, 'rb')
		msds = rd.MSDataSet()
		msds.load_file(sample_file2, header=True)

		peptide_filename = "./Hs_test.pepDict"
		peptide_file = open(peptide_filename)

		dm_peptide = nu.normalize_peptide_count(msds.get_data_matrix(), msds.get_id_dict(), peptide_file)
		
		print dm_peptide
		np.testing.assert_almost_equal( dm_peptide[4,3], 0.450261780104712) #kdrew: 86.0/191
		np.testing.assert_almost_equal( dm_peptide[5,2], 0.03333333333333333) #kdrew: 2.0/60


if __name__ == "__main__":
	unittest.main()


