#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.hierarchical_clustering as hc
import protein_complex_maps.physical_interaction_clustering as pc
import protein_complex_maps.cluster_util as cu

import numpy as np


class ClusterTest(unittest.TestCase):

	def setUp(self,):

		sample_filename = "./test_data.txt"
		sample_file = open(sample_filename, 'rb')
		self.data_matrix, self.name_list, self.fraction_list = rd.read_datafile(sample_file)

		#print self.data_matrix

		sample_filename2 = "./test_data2.txt"
		sample_file2 = open(sample_filename2, 'rb')
		self.data_matrix2, self.name_list2, self.fraction_list2 = rd.read_datafile(sample_file2)

		self.sample_file = open(sample_filename, 'rb')
		self.sample_file2 = open(sample_filename2, 'rb')
		#print self.data_matrix2

		self.msds = rd.MSDataSet()
		self.msds.load_file(self.sample_file, header=True)
		self.msds.load_file(self.sample_file2, header=True)

		self.msds.map_ids("ENSEMBL_ID", "ACC")

	def testHcluster(self,):

		fulldm = self.msds.get_data_matrix()
		print "get_data_matrix"
		print fulldm

		Y, Y2, D = hc.runCluster(fulldm)
		print Y
		assert (Y[0][0] == 4)
		assert (Y[0][1] == 7)
		assert (Y[1][0] == 8)
		assert (Y[1][1] == 9)

	def testGetCluster(self,):

		print "testGetCluster"
		fulldm = self.msds.get_data_matrix()
		print fulldm

		Y, Y2, D = hc.runCluster(fulldm)
		print Y
		print len(Y)
		cids = cu.get_cluster(Y, 13 )
		print cids
		assert (cids[0] == 4.0)
		assert (cids[1] == 7.0)
		cids = cu.get_cluster( Y, 0 )
		print cids
		assert (cids[0] == 0.0)

		cids = cu.get_cluster( Y, 24 )
		print cids
		assert ( len(cids) == 13 )

		try: 
			cids = cu.get_cluster( Y, 25 )
		except cu.ClusterOutOfBounds:
			pass

	def testPhysicalInteractions(self,):

		fulldm = self.msds.get_data_matrix()
		Y, Y2, D = hc.runCluster(fulldm)

		#pc.interactions_over_random( self.msds.get_name2index(), Y ) 
 		acc_mapping = self.msds.get_mapping("ACC_list")
		print acc_mapping
		print len(acc_mapping)
		pc.interactions_over_random( acc_mapping, Y ) 


if __name__ == "__main__":
	unittest.main()


