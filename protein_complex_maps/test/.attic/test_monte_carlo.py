#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.monte_carlo as mc
import protein_complex_maps.bicluster.bicluster as bc
import numpy as np

def score_func(matrix):
	return -1.0*matrix.shape[0]

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

		assert(self.montecarlo.iterations() == 1)
		assert(self.montecarlo.accept_rate() == 1.0)
		assert(self.montecarlo.current_score() == -2.0)

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		assert(self.montecarlo.iterations() == 2)
		assert(self.montecarlo.accept_rate() == 0.5)
		assert(self.montecarlo.current_score() == -2.0)


		self.bicluster1.add_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		assert(self.montecarlo.iterations() == 3)
		assert(self.montecarlo.accept_rate() == 2.0/3)
		assert(self.montecarlo.current_score() == -3.0)

		self.bicluster1.remove_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print "bclust: %s" % bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		assert(self.montecarlo.iterations() == 4)
		assert(self.montecarlo.accept_rate() == 2.0/4)
		assert(self.montecarlo.current_score() == -3.0)

		self.montecarlo.temp(1000.0)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print bclust.rows()
		print self.montecarlo.result_history()
		print self.montecarlo.score_history()
		print self.montecarlo.score_diff_history()
		print ""

		assert(self.montecarlo.iterations() == 5)
		assert(self.montecarlo.accept_rate() == 3.0/5)
		assert(self.montecarlo.current_score() == -2.0)

		#assert( score == 18008.325 ) #kdrew: divided by number of columns and rows


if __name__ == "__main__":
	unittest.main()


