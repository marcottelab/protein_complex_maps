
import protein_complex_maps.plots.plot_bicluster as pb
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.score_util as su
import protein_complex_maps.monte_carlo as mc
import protein_complex_maps.random_sampling_util as rsu

#kdrew: this class is used for seeding, generating, and evaluating biclusters
class BiclusterGenerator(object):

	def __init__( self, iterations=1000, random_module=None ): 

		if random_module == None:
			try:
				import numpy.random as random 
			except ImportError:
				self.random_module = None
			else:
				self.random_module = random
		else:
			self.random_module = random_module


		self.starting_temperature = 0.00000001
		self.montecarlo = mc.MonteCarlo(multiple_dot, temp=self.starting_temperature, random_module=random_module)
		self.iterations = iterations
		self.biclusters = []

	def random_seed(self, data_matrix):

		total_rows, total_columns = data_matrix.shape

		random_number_of_rows = self.random_module.random_integers(0,data_matrix.shape[0]-1)
		random_number_of_columns = self.random_module.random_integers(0,data_matrix.shape[1]-1)

		random_rows = self.random_module.choice(xrange(0,total_rows), random_number_of_rows, replace=False)
		random_columns = self.random_module.choice(xrange(0,total_columns), random_number_of_columns, replace=False)

		random_seeded_bicluster = bc.Bicluster(rows = random_rows, cols = random_columns, random_module=self.random_module)

		return random_seeded_bicluster

		#kdrew: temporary seed bicluster
		#bicluster1 = bc.Bicluster(rows = [3,9], cols = [20,22], random_module=self.random_module)
		#bicluster1 = bc.Bicluster(rows = [3,10,13,14,22,26], cols = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70], random_module=self.random_module)


	def generator(self,data_matrix, seed_bicluster=None):

		if seed_bicluster == None:
			#kdrew: if not seeded, use random
			bicluster1 = self.random_seed(data_matrix)
		else:
			#kdrew: copy seed bicluster to working bicluster
			bicluster1 = bc.Bicluster(rows=seed_bicluster.rows(), cols=seed_bicluster.columns(), random_module=self.random_module)
		

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
				print "row %s in bicluster" % (random_row,)
				bicluster1.remove_row(random_row)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			else:
				print "row %s out of bicluster" % (random_row,)
				bicluster1.add_row(random_row)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			if random_column in bicluster1.columns():
				print "column %s in bicluster" % (random_column,)
				bicluster1.remove_column(random_column)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)

			else:
				print "column %s out of bicluster" % (random_column,)
				bicluster1.add_column(random_column)
				bicluster1 = self.montecarlo.boltzmann(data_matrix, bicluster1)


			print self.montecarlo.result_history()
			print self.montecarlo.score_history()
			print self.montecarlo.score_diff_history()

		self.biclusters.append(self.montecarlo.lowscore_bicluster())


		self.evaluate(data_matrix, len(self.biclusters)-1)
		
		return self.montecarlo.lowscore_bicluster()


	#kdrew: calculates # of std away given bicluster score is from mean of randomly sampled biclusters 
	def evaluate(self, data_matrix, bc_index, plot=False): 

		score = multiple_dot( data_matrix, self.biclusters[bc_index] )

		lcb_rows = self.biclusters[bc_index].rows()
		lcb_columns = self.biclusters[bc_index].columns()
		dist_dict = rsu.random_sampling_score_distribution( data_matrix, su.multiple_dot, rows=lcb_rows, columns=lcb_columns, sample_module = self.random_module )
		for t in dist_dict.keys():
			zscore = abs(dist_dict[t].mean() - score)/dist_dict[t].std()
			print "random %s, mean: %s, std: %s, zscore: %s" % (t, dist_dict[t].mean(), dist_dict[t].std(), zscore)

		if plot:
			pb.plot_score_distribution(dist_dict, score=score, savefilename="/home/kdrew/public_html/test/bicluster_randomscore_plot.pdf")

def multiple_dot(data_matrix, bicluster):
	return -1.0*su.multiple_dot(bicluster.get_submatrix(data_matrix))





