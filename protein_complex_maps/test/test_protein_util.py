#!/usr/bin/python

# Unit tests for protein util 

import unittest
import MySQLdb
import protein_complex_maps.util.protein_util as pu
import numpy as np


class ProteinUtilTest(unittest.TestCase):

    def setUp(self,):
        self.protein_id = ["Q6N089","Q53HL2"]

    def testLength(self,):
        length = pu.get_length_uniprot(self.protein_id)
        assert( length["Q6N089"] == 472 )
        assert( length["Q53HL2"] == 280 )

    def testOrtholog(self,):
		
        try:
            ortholog_map = pu.get_ortholog(["Q6P2Q9","Q92736"], "Hsapiens", "Celegans")
        except MySQLdb.OperationalError as e:
            print "Test requires database: Is database set up?"
            print e
            return

        print ortholog_map
        assert( ortholog_map["Q6P2Q9"] == "P34369")
        assert( ortholog_map["Q92736"] == "I2HAA6")

		#ortholog_map2 = pu.get_ortholog(["ENSG00000197102","ENSG00000174231"], "Hsapiens", "Celegans", version="", database="inparanoid_blake")
		#print ortholog_map2
		#assert( ortholog_map2["ENSG00000197102"] == "T21E12.4")
		#assert( ortholog_map2["ENSG00000174231"] == "C50C3.6")


if __name__ == "__main__":
	unittest.main()


