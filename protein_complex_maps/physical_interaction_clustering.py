
import logging
import numpy as np
import os.path
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy
import pylab
import MySQLdb
import itertools as it
from scipy.stats.stats import pearsonr

import protein_complex_maps.normalization_util as nu
import protein_complex_maps.read_data as rd
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.plots.plot_profile as pp
import protein_complex_maps.bicluster_generator as bg
import protein_complex_maps.annealer as anl

import protein_complex_maps.protein_util as pu



import argparse
import pickle
import scipy.cluster.hierarchy as sch


logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

	parser = argparse.ArgumentParser(description="Physical interaction clusters fractionation data")
	#parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
	#					help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")

	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--complexe", action="store", dest="complexe", nargs='+', required=False, 
						help="Complex name from sql database in which to analyze")
	parser.add_argument("--input_complex_list", action="store", dest="complex_filename", required=False, 
						help="Filename of complex protein identifiers")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")

	args = parser.parse_args()

	#msds = pickle.load( open( args.msds_filename, "rb" ) )

	#if args.proteins != None:
	#	data_set, new_id_map = msds.get_subdata_matrix(args.proteins) 

	#elif args.complexe != None:
	#	db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'stoichiometry')
	#	cursor = db.cursor()
	#	try:
	#		cursor.execute("select distinct p.proteinid from complex as c, complex_association as ca, protein as p where p.id = ca.protein_key and ca.complex_key = c.id and c.name = '%s'" % (name) )
	#		complex_protein_ids = cursor.fetchall()
	#		complex_protein_ids = [ele for tupl in complex_protein_ids for ele in tupl]
	#		data_set, new_id_map = msds.get_subdata_matrix(args.proteins) 


	#	except MySQLdb.Error, e:
	#		print e

	#else:
	#	data_set = msds.get_data_matrix()
	#	new_id_map = msds.get_name2index()



	D = create_interaction_matrix( args.proteins )
	Y = runCluster( D )
	plotDendrogram(Y, D, args.plot_filename, args.proteins)
	#clusters = evaluateClusters(msds, data_set, new_id_map, fit_model, args.plot_filename, plot_threshold=args.plot_threshold)

def interactions_over_random( protein_id_map, Y ):
	print "in interactions_over_random"
	print Y
	print len(Y)
	print Y[-1]
	highest_cluster_id = len(Y) + len(protein_id_map) - 1
	print "highest_cluster_id: %s" % highest_cluster_id
	print "numOfProteins: %s" % len(protein_id_map)

	biogrid_map = pu.map_protein_ids( list(it.chain(*protein_id_map.values()))  , "ACC+ID", 'BIOGRID_ID' )

	#kdrew: clean up protein_id_map so only matches to biogrid_map are present
	for i in protein_id_map:
		not_found_flag = True
		for pid in protein_id_map[i]:
			#print "biogrid mapping: %s : %s" % (pid, biogrid_map[pid],)
			if len(biogrid_map[pid]) > 0:
				protein_id_map[i] = [pid]
				not_found_flag = False

		if not_found_flag:
			protein_id_map[i] = [pid]


	ior = IntOverRand(protein_id_map, Y, biogrid_map)

	#interactions_over_random_recursive( protein_id_map, Y, highest_cluster_id, biogrid_map )
	ior.interactions_over_random_recursive( highest_cluster_id )

	#print ior.zscores
	#print ior.distances

	return ior.zscores, ior.distances

class IntOverRand(object):
	def __init__(self, protein_id_map, Y, biogrid_map, nRandSamples=1000, sampling_type="all", bg_cutoff=0.5):
		self.protein_id_map = protein_id_map
		self.Y = Y
		self.biogrid_map = biogrid_map
		self.nRandSamples = nRandSamples
		self.sampling_type = sampling_type
		#kdrew: percent of background from which to sample in order to be valid
		self.bg_cutoff = bg_cutoff

		self.zscores = dict()
		self.distances = dict()

		protein_ids = list(it.chain(*self.protein_id_map.values() ))
		self.Dtotal = create_interaction_matrix(protein_ids , self.biogrid_map )
		#kdrew: add a bit of noise
		self.Dtotal = nu.sample_noise(self.Dtotal, sample_module=np.random.normal, scale=0.001)
		print "Dtotal: %s" % self.Dtotal

		self.missing_from_biogrid = 0
		for pid in protein_ids:
			if len(biogrid_map[pid]) == 0:
				self.missing_from_biogrid += 1

		print "missing from biogrid: %s out of %s" % (self.missing_from_biogrid, len(protein_ids))



	#def interactions_over_random_recursive( protein_id_map, Y, cluster_id, biogrid_map, nRandSamples=1000 ):
	def interactions_over_random_recursive( self, cluster_id ):

		print "in interactions_over_random recursive, %s" % cluster_id

		#kdrew: leaf node
		if cluster_id < len(self.protein_id_map):
			return [self.protein_id_map[cluster_id]]

		Y_id = cluster_id - len(self.protein_id_map)

		#kdrew: get all protids in subcluster i and j
		#kdrew: left subcluster
		proteins_lcluster = self.interactions_over_random_recursive( self.Y[Y_id][0] )
		
		#kdrew: right subcluster
		proteins_rcluster = self.interactions_over_random_recursive( self.Y[Y_id][1] )

		print "cluster_id %s Y_id %s" % (cluster_id, Y_id)
		print proteins_lcluster
		print proteins_rcluster

		#kdrew: unpack list
		protein_ids = list(it.chain(*proteins_lcluster+proteins_rcluster))

		#kdrew: get interactions for protids 
		D = create_interaction_matrix( protein_ids, self.biogrid_map )
		print D

		c_boundary = len(proteins_lcluster)
		cr_boundary = len(proteins_rcluster)
		lcluster_sum = D[:c_boundary+1,:c_boundary+1].sum()
		rcluster_sum = D[c_boundary:,c_boundary:].sum()

		rand_lcluster_sum_list = []
		rand_rcluster_sum_list = []
		#kdrew: for i in nRandSamples:
		for i in range(self.nRandSamples):
			#kdrew: draw random sample of protids of size len(protids[i]) from protids[i] + protids[j]

			#kdrew: shuffle is not what we want here, want to shuffle the ids in each cluster set
			#np.random.shuffle(D)

			#kdrew: samples from all proteins in whole complex
			if self.sampling_type == "all":
				rand_order = range(self.Dtotal.shape[0])
				np.random.shuffle(rand_order)
				#print "rand_order: %s" % rand_order
				Drand = self.Dtotal[:,rand_order][rand_order]
				rand_lcluster_sum = Drand[:c_boundary+1,:c_boundary+1].sum()
				rand_rcluster_sum = Drand[:cr_boundary+1,:cr_boundary+1].sum()

			#kdrew: samples from only the proteins in both the left and right clusters, all other proteins are not considered
			else:
				rand_order = range(D.shape[0])
				np.random.shuffle(rand_order)
				#print "rand_order: %s" % rand_order
				Drand = D[:,rand_order][rand_order]

				#print Drand

				#kdrew: get interactions for random protids 
				rand_lcluster_sum = Drand[:c_boundary+1,:c_boundary+1].sum()
				rand_rcluster_sum = Drand[c_boundary:,c_boundary:].sum()
				
			rand_lcluster_sum_list.append(rand_lcluster_sum)
			rand_rcluster_sum_list.append(rand_rcluster_sum)

			#print "rand lcluster_sum: %s" % lcluster_sum
			#print "rand rcluster_sum: %s" % rcluster_sum

		#kdrew: get indexes for left and right subclusters 
		l_cluster_id = self.Y[Y_id][0] 
		r_cluster_id = self.Y[Y_id][1] 
		lY_id = l_cluster_id - len(self.protein_id_map)
		rY_id = r_cluster_id - len(self.protein_id_map)

		l_cluster_zscore = (lcluster_sum - np.mean(rand_lcluster_sum_list)) / np.std(rand_lcluster_sum_list)
		r_cluster_zscore = (rcluster_sum - np.mean(rand_rcluster_sum_list)) / np.std(rand_rcluster_sum_list)

		print "lcluster_sum: %s zscore: %s" % (lcluster_sum, l_cluster_zscore)
		print "rcluster_sum: %s zscore: %s" % (rcluster_sum, r_cluster_zscore)

		c_boundary = len(proteins_lcluster)
		cr_boundary = len(proteins_rcluster)

		#kdrew: this makes sure the background from which we randomly sample is large enough
		l_bg_lrg_enough = len(proteins_lcluster) <= self.bg_cutoff * self.Dtotal.shape[0]
		r_bg_lrg_enough = len(proteins_rcluster) <= self.bg_cutoff * self.Dtotal.shape[0]

		#kdrew: if indices are valid, meaning they are not leaf nodes, 
		#kdrew: store distance metric and zscore
		if lY_id >= 0 and l_bg_lrg_enough: 
			ldist = self.Y[lY_id][2]
			self.distances[l_cluster_id] = ldist
			self.zscores[l_cluster_id] = l_cluster_zscore
		else:
			ldist = 0.0

		if rY_id >= 0 and r_bg_lrg_enough:
			rdist = self.Y[rY_id][2]
			self.distances[r_cluster_id] = rdist
			self.zscores[r_cluster_id] = r_cluster_zscore
		else:
			rdist = 0.0


		return proteins_lcluster + proteins_rcluster

def create_interaction_matrix( proteins, biogrid_map=None, map_from="ACC+ID", experimental_system=None):
	db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'biogrid')
	cursor = db.cursor()

	if type(experimental_system) is list:
		format_strings = ','.join(['"%s"'] * len(experimental_system))
		expsys_str = " and experimental_system in (%s)" % (format_strings,)
		expsys_str = expsys_str % (tuple(experimental_system))
	else:
		expsys_str = ""

	if biogrid_map == None:
		biogrid_map = pu.map_protein_ids( proteins, map_from, 'BIOGRID_ID' )

	D = scipy.zeros([len(proteins), len(proteins)])
	for i, prot in enumerate(proteins):
		for j, prot2 in enumerate(proteins):
			if len(biogrid_map[prot2]) == 0 or len(biogrid_map[prot]) == 0:
				#D[i,j] = -1
				continue
			print "prot: %s %s prot2: %s %s" % (prot, biogrid_map[prot], prot2, biogrid_map[prot2])
			cursor.execute("select * from Homo_sapiens_3_2_112 where throughput = 'Low Throughput' and BioGRID_ID_A = '%s' and BioGRID_ID_B = '%s' %s" % ( biogrid_map[prot][0], biogrid_map[prot2][0], expsys_str ) )
			rows = cursor.fetchall()
			print len(rows)
			D[i,j] = len(rows) != 0
			#D[i,j] = len(rows) 

			##kdrew: if i and j are same protein, mark them as being in the same physical interaction (this does not mean they are interacting with each other)
			#if i == j:
			#	D[i,j] = True
	
	return D


def runCluster(D): 

	Y = sch.linkage(D, method='centroid')

	return Y

#kdrew: code modified from stackoverflow
#http://stackoverflow.com/a/3011894
def plotDendrogram(Y, Y2, D, plot_filename, new_id_map):

	fig = pylab.figure(figsize=(8,8))
	ax1 = fig.add_axes([0.09, 0.1, 0.11, 0.6])
	dendrogram = sch.dendrogram(Y, orientation='right')
	print dendrogram['leaves']
	print [new_id_map[z] for z in dendrogram['leaves']]
	ax1.set_yticklabels([new_id_map[z] for z in dendrogram['leaves']])

	print new_id_map

	ax2 = fig.add_axes([0.3, 0.78, 0.6, 0.2])
	dendrogram2 = sch.dendrogram(Y2)
	ax2.set_xticklabels([new_id_map[z] for z in dendrogram2['leaves']], rotation='vertical', size='small')

	axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
	idx1 = dendrogram['leaves']
	idx2 = dendrogram2['leaves']
	D = D[idx1,:]
	D = D[:,idx2]
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
	axmatrix.set_xticks([])
	axmatrix.set_yticks([])


	axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])
	pylab.colorbar(im, cax=axcolor)
	

	print dendrogram

	#fig.show()
	fig.savefig(plot_filename)

	#print model.rows_
	#print model.columns_
	#print new_id_map

	return D



if __name__ == "__main__":
	main()


