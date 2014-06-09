
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

def create_interaction_matrix( proteins, map_from="ACC+ID"):
	db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'biogrid')
	cursor = db.cursor()

	biogrid_map = pu.map_protein_ids( proteins, map_from, 'BIOGRID_ID' )

	D = scipy.zeros([len(proteins), len(proteins)])
	for i, prot in enumerate(proteins):
		for j, prot2 in enumerate(proteins):
			print "prot: %s %s prot2: %s %s" % (prot, biogrid_map[prot], prot2, biogrid_map[prot2])
			cursor.execute("select * from Homo_sapiens_3_2_112 where throughput = 'Low Throughput' and BioGRID_ID_A = '%s' and BioGRID_ID_B = '%s'" % ( biogrid_map[prot][0], biogrid_map[prot2][0] ) )
			rows = cursor.fetchall()
			print len(rows)
			D[i,j] = len(rows) != 0

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


def evaluateClusters(msds, data_set, new_id_map, fit_model, plot_filename=None, plot_reordered=True, plot_threshold=0.0):

	bicluster_list = []

	rsscore_obj = rsu.RandomSamplingScore(data_set, su.multiple_dot_neg, sample_module=np.random)

	for j, row in enumerate(fit_model.rows_):
		print "Bicluster %s:" % (j,)

		bicluster_ids = []
		for i, val in enumerate(row):
			if val:
				try:
					#print "id %s, val %s, new_id_map %s"  % (i, val, new_id_map[i], )
					print "%s"  % (new_id_map[i], )
					bicluster_ids.append(new_id_map[i])
				except KeyError:
					continue


		#kdrew: matplotlib gets whiny about subplots when only one protein in plot, also make sure there are columns in matrix
		if len(bicluster_ids) > 1 and np.any(fit_model.columns_[j]):
			bicluster_data_set, bicluster_id_map = msds.get_subdata_matrix(bicluster_ids) 

			bicluster_data_set = bicluster_data_set[np.ix_(range(len(bicluster_data_set)), fit_model.columns_[j])]



			score = rsscore_obj.zscore_all_neg(bicluster_data_set)
			print "initial score: %s" % (score,)

			row_list = list(np.ix_(fit_model.rows_[j])[0]) 
			column_list = list(np.ix_(fit_model.columns_[j])[0])
			bicluster_list.append((bicluster_ids, score, row_list, column_list))


			#MC_OPT = False
			#if MC_OPT:
			#	bcgen = bg.BiclusterGenerator(rsscore_obj.zscore_all_neg, iterations=250, starting_temperature = 0.0, random_module=np.random)
			#	quench_annealer = anl.QuenchAnnealer( bcgen.get_montecarlo(), quench_iteration=250 )
			#	bcgen.set_annealer(quench_annealer)
            #
			#	for i in xrange(1,10):
			#		logger = logging.getLogger()
			#		logger.disabled = True
			#		bicluster1 = bcgen.generator(bicluster_data_set)
            #
			#		eval_dict = bcgen.evaluate( bicluster_data_set, len(bcgen.biclusters)-1 )
			#		logger.disabled = False
            #
			#		logging.info(bicluster1.rows())
			#		logging.info(bicluster1.columns())
            #
			#		for t in eval_dict.keys():
			#			logging.info("random %s, mean: %s, std: %s, zscore: %s" % ( t, eval_dict[t]['mean'], eval_dict[t]['std'], eval_dict[t]['zscore'] ))

			bc_msds = rd.MSDataSet()
			bc_msds.set_data_matrix(bicluster_data_set)
			bc_msds.set_id_dict({v:k for k, v in bicluster_id_map.items()})

			#print bc_msds.get_id_dict()

			if plot_filename != None and score < plot_threshold:
				plot_fname = os.path.splitext(plot_filename)[0]+'.'+str(j)+os.path.splitext(plot_filename)[1]
				#print bicluster_ids
				if plot_reordered:
					reordered_data_set, reordered_map = msds.reordered_data_matrix(row_list, column_list, data_set, new_id_map)
					numOfRows = sum(fit_model.rows_[j])
					numOfColumns = sum(fit_model.columns_[j])
					#kdrew: since we reordered, highlight the first # of rows and columns
					pp.plot_profile_dataset( reordered_data_set, reordered_map, savefilename=plot_fname, x_highlight= (0,numOfRows), y_highlight= (0,numOfColumns))
				else:
					pp.plot_profile( bc_msds, bicluster_ids, savefilename=plot_fname )


		print "\n"

	print new_id_map.values()
	if plot_filename != None:
		pp.plot_profile( msds, new_id_map.values(), total_occupancy = True, savefilename=plot_filename )

	return bicluster_list


if __name__ == "__main__":
	main()


