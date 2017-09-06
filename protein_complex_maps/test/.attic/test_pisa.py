#!/usr/bin/python

# Unit tests for pisa in pdb util 
import sys

import unittest
import MySQLdb
import protein_complex_maps.pdb_util as pdbu
import numpy as np
import Bio.PDB



class PISATest(unittest.TestCase):

    def setUp(self,):

        #kdrew: parse filename to get pdbid
        pisa_file = "./4cr2_pisa"
        pdbid = "4cr2"
        try:
            self.pisaInt = pdbu.PISA_Interfaces( pisa_file, pdbid=pdbid)
        except MySQLdb.OperationalError as e:
            print "Test requires database: Is database set up?"
            print e
            sys.exit()



    def testSurfaceArea(self,):
        surface_area = self.pisaInt.surface_area('1','2')
        print "surface_area: %s" % surface_area
        #np.testing.assert_almost_equal( surface_area, 773.6, decimal=1 )
        assert( surface_area == 773.6 )

        surface_area3 = self.pisaInt.surface_area('2','1')
        print "surface_area3: %s" % surface_area3
        #np.testing.assert_almost_equal( chain_Nterm_dist, 773.6, decimal=5 )
        assert( surface_area3 == 773.6 )

    def testSurfaceAreaByACC_Yeast(self,):
        surface_area5= self.pisaInt.surface_area_by_acc('P38624', 'P25043')
        print "surface_area5: %s" % surface_area5
        assert( surface_area5 == 773.6 )

        surface_area6 = self.pisaInt.surface_area_by_acc('P25043','P38624' )
        print "surface_area6: %s" % surface_area6
        #np.testing.assert_almost_equal( chain_Nterm_dist, 773.6, decimal=5 )

        assert( surface_area6 == 773.6 )

    def testSurfaceAreaByACC_Human(self,):
        surface_area2 = self.pisaInt.surface_area_by_acc('Q99436','P28072', base_species='Hsapiens')
        print "surface_area2: %s" % surface_area2
        assert( surface_area2 == 773.6 )

        surface_area4 = self.pisaInt.surface_area_by_acc('P28072','Q99436', base_species='Hsapiens')
        print "surface_area4: %s" % surface_area4
        assert( surface_area4 == 773.6 )


if __name__ == "__main__":
	unittest.main()


