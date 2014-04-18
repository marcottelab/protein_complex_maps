#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.read_data as rd
import protein_complex_maps.peptide_util as ppu
import numpy as np


class PeptideTest(unittest.TestCase):

	def setUp(self,):

		sample_filename2 = "./test_data2.txt"
		sample_file2 = open(sample_filename2, 'rb')
		msds = rd.MSDataSet()
		msds.load_file(sample_file2, header=True)

		peptide_filename = "./Hs_test.pepDict"
		peptide_file = open(peptide_filename)

		self.peptide_dict = ppu.read_peptide_dict(peptide_file, msds.get_id_dict().keys())


		peptide_list_filename = "./Hs_test.peplist"
		peptide_list_file = open(peptide_list_filename)

		self.peptide_count_dict = ppu.read_peptide_list(peptide_list_file)



	def testPeptideDict(self, ):
		
		#print self.peptide_dict
		#print self.peptide_dict['ENSG00000087191'] 
		assert( self.peptide_dict['ENSG00000087191'][0] == 'AGSGLRQYYLSK' )

	def testPeptideList(self, ):
		#print self.peptide_count_dict

		#for peptide in self.peptide_count_dict:
			#print peptide
			#print self.peptide_count_dict[peptide]

		assert (self.peptide_count_dict['DYETATLSEIK']['PH090902_HS3NE_HCW_P2B03'] == 1.0)


	def testPeptides(self,):

		#for protid in self.peptide_dict:
		#	for peptide in self.peptide_dict[protid]:
		#		try:
		#			print self.peptide_count_dict[peptide], protid, peptide
		#		except KeyError:
		#			continue


		protein_count_dict = ppu.protein_counts_by_peptide(self.peptide_dict, self.peptide_count_dict)
		#for prot in protein_count_dict:
		#	for fraction in protein_count_dict[prot]:
				#print fraction
				#print protein_count_dict[prot][fraction]


		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['DYETATLSEIK'] == 1.0 )
		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['KDDPLTNLNTAFDVAEK'] == 1.0 )
		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['RDQALTEEHAR'] == 1.0 )

	def testMSDM(self,):

		protein_count_dict = ppu.protein_counts_by_peptide(self.peptide_dict, self.peptide_count_dict)
		msds2 = rd.MSDataSet()
		#print "calling create_by_peptide_counts:"
		msds2.create_by_peptide_counts( protein_count_dict )
		#print msds2.get_data_matrix()
		#print msds2.get_id_dict()
		#print msds2.get_fraction_dict()

		assert( msds2.get_data_matrix()[msds2.get_id_dict()['ENSG00000010438'],msds2.get_fraction_dict()['PH090902_HS3NE_HCW_P2B05']] == 2.0)

		

if __name__ == "__main__":
	unittest.main()


