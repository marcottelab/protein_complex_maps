
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.score_util as su
import protein_complex_maps.monte_carlo as mc

class BiclusterGenerator(object):

	def __init__( self, iterations=1000, random_module=None ): 

		self.montecarlo = mc.MonteCarlo(multiple_dot, temp=0.0, random_module=random_module)
		self.iterations = iterations

		if random_module == None:
			try:
				import numpy.random as random 
			except ImportError:
				self.random_module = None
			else:
				self.random_module = random
		else:
			self.random_module = random_module

	def generator(self,data_matrix):

		#kdrew: temporary seed bicluster
		bicluster1 = bc.Bicluster(rows = [3,9], cols = [20,22], random_module=self.random_module)

		numRows, numCols = data_matrix.shape
		print numRows, numCols

		#kdrew: probably should test for some convergence
		for i in xrange( 1, self.iterations ):
			print "iteration: %s" % (i,)

			#kdrew: randomly pick row or column (inside or outside of bicluster)
			random_row = self.random_module.random_integers(0, numRows-1)
			random_column = self.random_module.random_integers(0, numCols-1)

			print "bicluster matrix: %s" % (bicluster1.get_submatrix(data_matrix))

			if random_row in bicluster1.rows():
				bicluster1.remove_row(random_row)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			else:
				bicluster1.add_row(random_row)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			if random_column in bicluster1.columns():
				bicluster1.remove_column(random_column)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			else:
				bicluster1.add_column(random_column)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)


			print self.montecarlo.result_history()

		return self.montecarlo.lowscore_bicluster()

def multiple_dot(data_matrix, bicluster):
	return su.multiple_dot(bicluster.get_submatrix(data_matrix))





