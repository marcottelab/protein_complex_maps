#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import numpy as np


class ReadDataTest(unittest.TestCase):

	def setUp(self,):

		sample_filename = "./test_data.txt"
		sample_file = open(sample_filename, 'rb')
		self.data_matrix, self.name_list, self.fraction_list = rd.read_datafile(sample_file)

		print self.data_matrix

		sample_filename2 = "./test_data2.txt"
		sample_file2 = open(sample_filename2, 'rb')
		self.data_matrix2, self.name_list2, self.fraction_list2 = rd.read_datafile(sample_file2)

		self.sample_file = open(sample_filename, 'rb')
		self.sample_file2 = open(sample_filename2, 'rb')
		print self.data_matrix2

		self.sample_Ce_filename = "./test_Ce_data.txt"
		self.sample_Ce_file = open(self.sample_Ce_filename,'rb')

	def testRead(self,):
		
		assert(self.data_matrix[np.ix_([0],[0])] == 0.0)
		assert(self.data_matrix[np.ix_([2],[0])] == 10.0)
		assert(self.data_matrix[np.ix_([2],[1])] == 100.0)
		assert(self.data_matrix[np.ix_([11],[2])] == 20.0)
		
		assert(self.data_matrix2[np.ix_([0],[0])] == 4.0)
		assert(self.data_matrix2[np.ix_([2],[0])] == 0.0)
		assert(self.data_matrix2[np.ix_([2],[1])] == 0.0)
		assert(self.data_matrix2[np.ix_([2],[3])] == 121.0)
		assert(self.data_matrix2[np.ix_([10],[3])] == 5.0)
		#np.testing.assert_almost_equal( score, -1.73220495821 )
		#np.testing.assert_almost_equal( score, 0.000343732988528)

	def testConcat(self,):
		print "testConcat"
		c_dmat, c_nlist = rd.concat_data_matrix(self.data_matrix, self.name_list, self.data_matrix2, self.name_list2)

		print c_dmat
		print c_nlist

		index1 = c_nlist.index('ENSG00000072110')
		assert(c_dmat[np.ix_([index1],[0])] == 25.0)
		assert(c_dmat[np.ix_([index1],[7])] == 30.0)

		index2 = c_nlist.index('ENSG00000100519')
		assert(c_dmat[np.ix_([index2],[3])] == 105.0)
		assert(c_dmat[np.ix_([index2],[7])] == 0.0)

	def testCeMsDataSet(self,):
		msds = rd.MSDataSet()
		msds.load_file(self.sample_Ce_file, header=True)

		msds.map_ids_by_genename(organism="Caenorhabditis+elegans")

		mapped_ids = msds.get_id_dict()

		assert(mapped_ids["C18H9.2"] == 6)
		assert(mapped_ids["Q94242"] == 5)
		assert(mapped_ids["Q94241"] == 2)
		assert(mapped_ids["F55A4.7"] == 5)


	def testMsDataSet(self,):
		msds = rd.MSDataSet()
		msds.load_file(self.sample_file, header=True)
		msds.load_file(self.sample_file2, header=True)

		fulldm = msds.get_data_matrix()
		print "get_data_matrix"
		print fulldm

		#kdrew: removed this functionality
		#partdm = msds.get_data_matrix(names=["ENSG00000072110","ENSG00000075914","ENSG00000087191"])
		#print "get_data_matrix"
		#print partdm
		#assert(partdm.shape[0] == 3)

		msds.map_ids("ENSEMBL_ID", "ACC")
		mapped_ids = msds.get_id_dict()
		assert(mapped_ids["Q15024"] == 9)
		assert(mapped_ids["G3V380"] == 11)
		assert(mapped_ids["J3QRR3"] == 6)       
		assert(mapped_ids["ENSG00000072110"] == 11)

		frac_map = msds.get_fraction_dict()
		assert(frac_map['PH090807_HS3NE_HCW_P1A03'] == 2)
		assert(frac_map['PH090807_HS3NE_HCW_P1A03a'] == 8)

	def testReorder(self,):
		msds = rd.MSDataSet()
		msds.load_file(self.sample_file, header=True)
		msds.load_file(self.sample_file2, header=True)

		fulldm = msds.get_data_matrix()
		print "testing Reorder"
		print fulldm

		reordered_dm, new_map = msds.reordered_data_matrix([5,6,2], [7,8,3])

		print reordered_dm

		print "name2index: %s" % msds.get_name2index()
		print "new_map: %s" % new_map


		assert msds.get_name2index()[5] == new_map[0]
		assert msds.get_name2index()[6] == new_map[1]
		assert msds.get_name2index()[2] == new_map[2]
		assert msds.get_name2index()[0] == new_map[3]

		assert fulldm[5,7] == reordered_dm[0,0]
		assert fulldm[6,7] == reordered_dm[1,0]
		assert fulldm[6,8] == reordered_dm[1,1]
		assert fulldm[2,3] == reordered_dm[2,2]
		assert fulldm[0,0] == reordered_dm[3,3]


if __name__ == "__main__":
	unittest.main()


