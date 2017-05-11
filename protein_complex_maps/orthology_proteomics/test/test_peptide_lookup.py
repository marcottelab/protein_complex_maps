#!/usr/bin/python


# Unit tests for new clustering analysis functionality

import sys
sys.path.append('../')
sys.path.append('.')
import unittest
import scripts.trypsin as tr
import scripts.define_grouping as dg
import scripts.get_elution_profiles as gep

from Bio import SeqIO

#import protein_complex_maps.read_data as rd
#import protein_complex_maps.peptide_util as ppu
import numpy as np
import pandas as pd

class LookupTest(unittest.TestCase):

        def setUp(self,):
                sample_proteome_filename = "test/test_proteome.fasta"
                sample_proteome = list(SeqIO.parse(sample_proteome_filename, "fasta"))
 
                self.sample_orthology = "test/test_orthology.tab"
                #sample_orthology = open(sample_orthology_filename, "r")
                self.sample_peptides = "test/test_peptides.csv"


                self.testprotein = sample_proteome[0].seq 
    
                self.contam_file = "test/contam_benzo_peptides.csv"
                self.sample_elution = "test/test_elution.csv"

               

#        def testDigest(self,):
#                #Make sure the right test protein is being chopped
#                assert(str(self.testprotein) =="MALMMKSSASLKAVSAGRSRRAVVVRAGKYDEELIKTAGTVASKGRGILAMDESNATCGKRLDSIGVENTEENRRAYRELLVTAPGLGQYISGAILFEETLYQSTASGKKFVDVMKEQNIVPGIKVDKGLVPLSNTNGESWCMGLDGLDKRCAEYYKAGARFAKWRSVVSIPHGPSIIAARDCAYGLARYAAIAQNAGLVPIVEPEVLLDGEHDIDRCLEVQEAIWAETFKYMADNKVMFEGILLKPAMVTPGADCKNKAGPAKVAEYTLKMLRRRVPPAVPGIMFLSGGQSELESTLNLNAMNQSPNPWHVSFSYARALQNTVLKTWQGKPENVQAAQAALLKRAKANSDAQQGKYDATTEGKEAAQGMYEKGYVY")
#                #Check that the right number of peptides are being made
#
#                self.peptides = tr.TRYPSIN(self.testprotein, 2)
#
#                print(len(self.peptides))
#                assert(len(self.peptides)== 108)
#
#        def testUniqueProteins(self,):
#                peps = dg.format_peptides(self.sample_peptides, 'test')
#
#                prot_uniq_peps = dg.protein_uniq_peptides(peps, 'test', 'test/')
#              
#               #Should not be unique to any protein
#                testrow1 = prot_uniq_peps[prot_uniq_peps.Peptide== 'AAAAAAKAJEQQKAEJESAQQQJESAR']
#                print(testrow1)
#                assert(len(testrow1)==0)
#
#               #Should not be unique to a single protein
#                testrow2 = (prot_uniq_peps[prot_uniq_peps.Peptide== 'AQAESCRHVWR'])
#                assert(len(testrow2)==0)
#
#                #Should get one unique protein
#                testrow3 = (prot_uniq_peps[prot_uniq_peps.Peptide== 'AAAAHAJTGVJKQASGRK'])
#                assert(len(testrow3)==1)
#                assert(testrow3['ProteinID'].values[0] == 'tr|A8J4W6|A8J4W6_CHLRE') 
#        def testGrouping(self,):
#
#                peps = dg.format_peptides(self.sample_peptides, 'test')
#           
#                group_uniq_peps = dg.define_grouping(self.sample_orthology, peps, 'test', 'euNOG', 'test/')
#
#                print(group_uniq_peps)
#
#                #Should not be unique to any protein
#                testrow1 = (group_uniq_peps[group_uniq_peps.Peptide== 'AAAAAAKAJEQQKAEJESAQQQJESAR'])
#                assert(len(testrow1)==0)
#
#                #Should get one unique groupein
#                testrow2 = (group_uniq_peps[group_uniq_peps.Peptide== 'AQAESCRHVWR'])
#                assert(len(testrow2)==1)
#                assert(testrow2['ID'].values[0] == 'KOG0001') 
#                 
#
#                testrow3 = (group_uniq_peps[group_uniq_peps.Peptide== 'AAAAHAJTGVJKQASGRK'])
#                assert(len(testrow3)==1)
#                assert(testrow3['ID'].values[0] == 'KOG0027') 
 
        def testElution(self,):

                gep.create_tables('test', 'euNOG', 'testing', self.sample_elution, self.sample_peptides, self.contam_file)





 


#
#
#	def setUp(self,):
#
#		sample_filename2 = "./test_data2.txt"
#		sample_file2 = open(sample_filename2, 'rb')
#		msds = rd.MSDataSet()
#		msds.load_file(sample_file2, header=True)
#
#		peptide_filename = "./Hs_test.pepDict"
#		peptide_file = open(peptide_filename)
#
#                self.
#
#
#
#		self.peptide_dict = ppu.read_peptide_dict(peptide_file, msds.get_id_dict().keys())
#
#		peptide_file = open(peptide_filename)
#		self.peptide_dict2 = ppu.read_peptide_dict(peptide_file)
#		print "self.peptide_dict2:"
#		for key in self.peptide_dict2:
#			print key
#			for pep in self.peptide_dict2[key]:
#				print pep
#
#			print ''
#
#
#		peptide_list_filename = "./Hs_test.peplist"
#		peptide_list_file = open(peptide_list_filename)
#
#		self.peptide_count_dict = ppu.read_peptide_list(peptide_list_file)
#
#
#
#	def testPeptideDict(self, ):
#		
#		#print self.peptide_dict
#		#print self.peptide_dict['ENSG00000087191'] 
#		assert( self.peptide_dict['ENSG00000087191'][0] == 'AGSGLRQYYLSK' )
#
#	def testPeptideList(self, ):
#		#print self.peptide_count_dict
#
#		#for peptide in self.peptide_count_dict:
#			#print peptide
#			#print self.peptide_count_dict[peptide]
#
#		assert (self.peptide_count_dict['DYETATLSEIK']['PH090902_HS3NE_HCW_P2B03'] == 1.0)
#
#
#	def testPeptides(self,):
#
#		#for protid in self.peptide_dict:
#		#	for peptide in self.peptide_dict[protid]:
#		#		try:
#		#			print self.peptide_count_dict[peptide], protid, peptide
#		#		except KeyError:
#		#			continue
#
#
#		protein_count_dict = ppu.protein_counts_by_peptide(self.peptide_dict, self.peptide_count_dict)
#		#for prot in protein_count_dict:
#		#	for fraction in protein_count_dict[prot]:
#				#print fraction
#				#print protein_count_dict[prot][fraction]
#
#
#		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['DYETATLSEIK'] == 1.0 )
#		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['KDDPLTNLNTAFDVAEK'] == 1.0 )
#		assert( protein_count_dict['ENSG00000072110']['PH090902_HS3NE_HCW_P2B03']['RDQALTEEHAR'] == 1.0 )
#
#	def testMSDM(self,):
#
#		protein_count_dict = ppu.protein_counts_by_peptide(self.peptide_dict, self.peptide_count_dict)
#		msds2 = rd.MSDataSet()
#		#print "calling create_by_peptide_counts:"
#		msds2.create_by_peptide_counts( protein_count_dict )
#		#print msds2.get_data_matrix()
#		#print msds2.get_id_dict()
#		#print msds2.get_fraction_dict()
#
#		dmat = msds2.get_data_matrix()
#		print "dmat sum : %s" % dmat.sum()
#		print "\n"
#
#		assert( msds2.get_data_matrix()[msds2.get_id_dict()['ENSG00000010438'],msds2.get_fraction_dict()['PH090902_HS3NE_HCW_P2B05']] == 2.0)
#
#	def testMSDM_threshold(self,):
#		protein_count_dict = ppu.protein_counts_by_peptide(self.peptide_dict, self.peptide_count_dict)
#		msds2 = rd.MSDataSet()
#		msds2.create_by_peptide_counts( protein_count_dict, threshold=2 )
#
#		dmat = msds2.get_data_matrix()
#		print "dmat sum threshold: %s" % dmat.sum()
#		print "\n"
#		

if __name__ == "__main__":
	unittest.main()


