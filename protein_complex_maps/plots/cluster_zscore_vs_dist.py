
import logging
import numpy as np
import os.path
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy
import pylab
import itertools as it
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.read_data as rd
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.plots.plot_profile as pp
import protein_complex_maps.bicluster_generator as bg
import protein_complex_maps.annealer as anl
import protein_complex_maps.physical_interaction_clustering as pic
import protein_complex_maps.hierarchical_clustering as hc


import argparse
import pickle
import scipy.cluster.hierarchy as sch

import protein_complex_maps.protein_util as pu

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

	parser = argparse.ArgumentParser(description="Plot zscore vs distance from physical interactions in subcomplexes found by hierarchical clustering")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--complexe", action="store", dest="complexe", nargs='+', required=False, 
						help="Complex name from sql database in which to analyze")
	parser.add_argument("--input_complex_list", action="store", dest="complex_filename", required=False, 
						help="Filename of complex protein identifiers")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
						help="Ignore missing protein ids in msds")
	parser.add_argument("--cluster_method", action="store", dest="cluster_method", required=False, default="single",
						help="""Type of linkage clustering, 
						types found: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.cluster.hierarchy.linkage.html, 
						default: single""")

	#parser.add_argument("-j", action="store", dest="numOfProcs", required=False, default=1,
	#    				help="Number of processors to use, default=1")

	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )

	if args.proteins != None:
		data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 

	elif args.complexe != None:
		db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'stoichiometry')
		cursor = db.cursor()
		try:
			cursor.execute("select distinct p.proteinid from complex as c, complex_association as ca, protein as p where p.id = ca.protein_key and ca.complex_key = c.id and c.name = '%s'" % (name) )
			complex_protein_ids = cursor.fetchall()
			complex_protein_ids = [ele for tupl in complex_protein_ids for ele in tupl]
			data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 


		except MySQLdb.Error, e:
			print e

	else:
		data_set = msds.get_data_matrix()
		new_id_map = msds.get_name2index()

	zscore_list, distance_list = compute_cluster_zscores( args.proteins, data_set, new_id_map, msds, args.cluster_method )
	plot_cluster_zscore_vs_dist( zscore_list, distance_list, args.plot_filename )

def compute_cluster_zscores( proteins, data_set, new_id_map, msds, cluster_method ):

	Y, Y2, D = hc.runCluster( data_set, cluster_method=cluster_method )

	acc_mapping = msds.get_mapping("ACC_list")
	#print acc_mapping
	#print len(acc_mapping)

	map_ids = msds.get_id_dict()
	#print map_ids
	#print new_id_map
	new_map_ids = {v:k for k, v in new_id_map.items()}
	#print new_map_ids
	new_acc_mapping = dict()
	for protid in proteins:
		i = map_ids[protid]
		#print i
		acc_list = acc_mapping[i]
		#print acc_list
		new_i = new_map_ids[protid]
		#print new_i
		new_acc_mapping[new_i] = acc_list

	#print new_acc_mapping 

	zscores, distances = pic.interactions_over_random( new_acc_mapping, Y ) 
	print zscores
	print distances

	zscore_list = []
	distance_list = []

	for k in zscores:
		zscore_list.append(zscores[k])
		distance_list.append(distances[k])

	return zscore_list, distance_list


def plot_cluster_zscore_vs_dist( zscore_list, distance_list, plot_filename ):

	fig = pylab.figure()
	pylab.scatter(distance_list, zscore_list)

	fig.savefig(plot_filename)


if __name__ == "__main__":
	main()


