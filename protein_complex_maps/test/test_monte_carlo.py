#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.monte_carlo as mc
import protein_complex_maps.bicluster.bicluster as bc
import numpy as np

def score_func(data_matrix, bicluster):
	return -1.0*len(bicluster.rows())

class MonteCarloTest(unittest.TestCase):

	def setUp(self,):

		np.random.seed(123)

		self.data_matrix = np.arange(40).reshape(4,10)

		self.bicluster1 = bc.Bicluster( rows = [1,2], cols = [4,5] ) 
		self.montecarlo = mc.MonteCarlo(score_func, temp=0.0, random_module=np.random)

	def testAcceptReject(self, ):

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		self.bicluster1.add_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		self.bicluster1.remove_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		self.montecarlo.temp(1000.0)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		#assert( score == 18008.325 ) #kdrew: divided by number of columns and rows


if __name__ == "__main__":
	unittest.main()


