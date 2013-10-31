#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.monte_carlo as mc
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.annealer as anl
import numpy as np

#kdrew: dummy score function that returns the size of bicluster
def score_func(matrix):
	return -1.0*matrix.shape[0]

class AnnealerTest(unittest.TestCase):

	def setUp(self,):

		np.random.seed(123)

		self.data_matrix = np.arange(40).reshape(4,10)

		self.bicluster1 = bc.Bicluster( rows = [1,2], cols = [4,5] ) 
		self.montecarlo = mc.MonteCarlo(score_func, temp=10.0, random_module=np.random)

	def testAnnealer(self, ):
		annealer =  anl.Annealer(self.montecarlo)

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		annealer.anneal()
		assert(self.montecarlo.temp() == 5.0)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		annealer.anneal()
		assert(self.montecarlo.temp() == 2.5)

	def testQuenchAnnealer(self, ):
		annealer =  anl.QuenchAnnealer( self.montecarlo, 2 )

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		annealer.anneal()
		assert(self.montecarlo.temp() == 10.0)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		annealer.anneal()
		assert(self.montecarlo.temp() == 0.0)
 
	def testRateAnnealer(self, ):
		#kdrew: make temp accept 75%
		annealer =  anl.RateAnnealer( self.montecarlo, rate=0.75 )

		self.montecarlo.temp(10.0)

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print self.montecarlo.temp()
		print self.montecarlo.accept_rate()
		#kdrew: accept causes annealer to reduce temp
		annealer.anneal()
		print self.montecarlo.temp()
		assert(self.montecarlo.temp() == 5.0)

		self.bicluster1.remove_row(1)
		self.bicluster1.remove_row(2)
		self.bicluster1.remove_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		#kdrew: rejection causes annealer to raise temp
		annealer.anneal()
		assert(self.montecarlo.temp() == 10.0)

	def testRampAnnealer(self, ):
		annealer =  anl.RampAnnealer( self.montecarlo, ramp=2 )

		self.montecarlo.temp(5.0)

		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print self.montecarlo.accept_rate()
		annealer.anneal()
		assert(self.montecarlo.temp() == 5.0)

		self.bicluster1.remove_row(1)
		self.bicluster1.remove_row(2)
		self.bicluster1.remove_row(3)
		bclust = self.montecarlo.boltzmann(self.data_matrix, self.bicluster1)
		print self.montecarlo.accept_rate()
		#kdrew: after two iterations the temp is set to original starting temp (see setup function)
		annealer.anneal()
		assert(self.montecarlo.temp() == 10.0)

if __name__ == "__main__":
	unittest.main()

 
