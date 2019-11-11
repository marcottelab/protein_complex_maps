#!/usr/bin/python

# Unit tests for protein util 

import unittest
import protein_complex_maps.util.protein_util as pu
import numpy as np


class ProteinUtilTest(unittest.TestCase):

    def setUp(self,):
        self.protein_id = ["Q6N089","Q53HL2"]

    def testLength(self,):
        length = pu.get_length_uniprot(self.protein_id)
        assert( length["Q6N089"] == 472 )
        assert( length["Q53HL2"] == 280 )

if __name__ == "__main__":
	unittest.main()


