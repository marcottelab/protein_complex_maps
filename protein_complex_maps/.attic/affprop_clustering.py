
import logging
import numpy as np
import os.path
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy
import scipy.sparse as sparse
import pylab
import itertools as it
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from sklearn.cluster import AffinityPropagation
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.read_data as rd
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.plots.plot_profile as pp
import protein_complex_maps.bicluster_generator as bg
import protein_complex_maps.annealer as anl
import protein_complex_maps.physical_interaction_clustering as pic
import protein_complex_maps.external.npeet.entropy_estimators as ee


import argparse
import pickle
import scipy.cluster.hierarchy as sch

import protein_complex_maps.protein_util as pu

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

	parser = argparse.ArgumentParser(description="Affinity propagation clusters fractionation data")
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
	parser.add_argument("--physical_plot_filename", action="store", dest="physical_plot_filename", required=False, default=None,
						help="Filename of physical interaction clustered plot")
	parser.add_argument("--plot_profile", action="store_true", dest="plot_profile", required=False, default=False,
						help="Plot profile instead of correlation matrix")
	parser.add_argument("--plot_just_dendrogram", action="store_true", dest="plot_just_dendrogram", required=False, default=False,
						help="Plot just dendrogram without correlation matrix")
	parser.add_argument("--pickle_filename", action="store", dest="pickle_filename", required=False, default=None,
						help="Filename of linkage object pickle, if set will pickle two linkages and correlation matrix in tuple")
	parser.add_argument("--sample_method", action="store", dest="sample_method", required=False, default=None,
						help="Sampling method to add noise to correlation calculations (poisson or normal)")
	parser.add_argument("--average_cnt", action="store", type=int, dest="average_cnt", required=False, default=0,
						help="Number of samples to compute average correlation")
	parser.add_argument("--total_occupancy", action="store_true", dest="total_occupancy", required=False, default=False,
						help="Flag to only plot columns where atleast one protein has greater than zero counts")
	parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
						help="Ignore missing protein ids in msds")
	parser.add_argument("--genenames", action="store_true", dest="genenames", required=False, default=False,
						help="Set labels to genenames")
	parser.add_argument("--cluster_method", action="store", dest="cluster_method", required=False, default="single",
						help="""Type of linkage clustering, 
						types found: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.cluster.hierarchy.linkage.html, 
						default: single""")
	parser.add_argument("--dependence_metric", action="store", dest="dependence_metric", required=False, default="pearson",
						help="""pearson or mutual_information
						default: pearson""")

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

	if args.genenames:
		genename_map = pu.get_genenames_uniprot( new_id_map.values() )
		print new_id_map
		print genename_map
		gene_id_map = dict()
		for i in xrange(len(data_set)):
			print i
			print genename_map[new_id_map[i]]
			#kdrew: sometimes no genename is returned for certain ids, default to original id
			if genename_map[new_id_map[i]] == None:
				gene_id_map[i] = new_id_map[i]
			else:
				gene_id_map[i] = genename_map[new_id_map[i]]

		new_id_map = gene_id_map

	if args.sample_method == "poisson":
		sample_module = np.random.poisson
	elif args.sample_method == "normal":
		sample_module = np.random.normal 
	else:
		sample_module = None

	af = runCluster( data_set, new_id_map )

	labels = af.labels_

	for i in range(labels.max() + 1):
		indices = np.where(labels == i)[0]
		print "Cluster %i: %s" % ((i+1), ','.join([ new_id_map[j] for j in indices]))

	print labels
	

#kdrew: both average_cnt and sample_module need to be set for average correlation with sample noise to be computed
#kdrew: consider moving correlation to another module
def runCluster(data_set, id_map, threshold=0.1): 

	data_set = np.nan_to_num(data_set)

	corrcoefMat = np.nan_to_num(np.corrcoef(data_set))
	for i in xrange(corrcoefMat.shape[0]):
		for j in xrange(corrcoefMat.shape[1]):
			print "%s\t%s\t%s" % (id_map[i], id_map[j], corrcoefMat[i,j],)

	#low_indices = corrcoefMat < threshold
	#corrcoefMat[low_indices] = 0.0
	#corrcoefMat_sparse = sparse.csc_matrix( corrcoefMat )

	#af = AffinityPropagation(affinity="precomputed").fit(corrcoefMat_sparse)
	
	return #af


if __name__ == "__main__":
	main()


