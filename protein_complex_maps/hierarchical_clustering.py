
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

import numpy.linalg as linalg

from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

import sklearn.metrics as skm

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

	parser = argparse.ArgumentParser(description="Hierarchical clusters fractionation data")
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
						help="""pearson or mutual_information or mutual_information_score
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
		print "new_id_map.values: %s" % (new_id_map.values(),)
		genename_map = pu.get_genenames_uniprot( new_id_map.values() )
		print new_id_map
		print genename_map
		gene_id_map = dict()
		for i in xrange(len(data_set)):
			print i
			#print genename_map[new_id_map[i]]
			#kdrew: sometimes no genename is returned for certain ids, default to original id
			try:
				if genename_map[new_id_map[i]] == None:
					gene_id_map[i] = new_id_map[i]
				else:
					gene_id_map[i] = genename_map[new_id_map[i]]
			except KeyError:
				gene_id_map[i] = new_id_map[i]

		new_id_map = gene_id_map

	if args.sample_method == "poisson":
		sample_module = np.random.poisson
	elif args.sample_method == "normal":
		sample_module = np.random.normal 
	else:
		sample_module = None

	Y, Y2, D = runCluster( data_set, args.average_cnt, sample_module, cluster_method=args.cluster_method, dependence_metric=args.dependence_metric )

        #Dinv = linalg.inv(D)
        #print Dinv

	
	if args.pickle_filename != None:
		pickle.dump((Y,Y2,D,new_id_map), open(args.pickle_filename, "wb"))

	if args.plot_profile:
		D = plotDendrogramProfile(data_set, Y, args.plot_filename, new_id_map, total_occupancy=args.total_occupancy, title=args.msds_filename)
	elif args.plot_just_dendrogram:
		D = plotJustDendrogram(Y, args.plot_filename, new_id_map, title=args.msds_filename)
	else:
		#kdrew: reordered correlation matrix gets returned
		D = plotDendrogram(Y, Y2, D, args.plot_filename, new_id_map, title=args.msds_filename)


	if args.physical_plot_filename != None:
		#kdrew: plot physical interactions
		pD = pic.create_interaction_matrix( args.proteins )
		pD = pic.plotDendrogram(Y, Y2, pD, args.physical_plot_filename, new_id_map)

		D_flat = np.array(D.reshape(-1)).flatten()
		pD_flat = np.array(pD.reshape(-1)).flatten()
		corr, pval = pearsonr(D_flat, pD_flat)
		print "pearsonr: %s pval: %s" % (corr, pval,)
		corr, pval = spearmanr(D_flat, pD_flat)
		print "spearmanr: %s pval: %s" % (corr, pval,)
		D_list = [[x] for x in D_flat.tolist()]
		pD_list = [[x] for x in pD_flat.tolist()]
		minfo= ee.mi(D_list, pD_list)
		print "mutual information: %s" % (minfo, )

#kdrew: both average_cnt and sample_module need to be set for average correlation with sample noise to be computed
#kdrew: consider moving correlation to another module
def runCluster(data_set, average_cnt=0, sample_module=None, scale=None, cluster_method="single", dependence_metric="pearson"): 

	data_set = np.nan_to_num(data_set)
	#data_set = nu.add_noise_over_columns(data_set)
	#data_set = nu.add_noise(data_set, 0.00000001)


	#kdrew: create distance matrix the size of proteins v proteins
	D = scipy.zeros([data_set.shape[0], data_set.shape[0]])

	for i in xrange(data_set.shape[0]):
		for j in xrange(data_set.shape[0]):
			if dependence_metric == "pearson":
				if average_cnt > 0 and sample_module != None:
					corr_sum = 0.0
					for k in range(average_cnt):
						a1 = np.nan_to_num(np.array(data_set[i,])[0])
						a2 = np.nan_to_num(np.array(data_set[j,])[0])

						if scale != None:
							a1 = nu.sample_noise(a1, sample_module, scale=scale)
							a2 = nu.sample_noise(a2, sample_module, scale=scale)
						else:
							a1 = nu.sample_noise(a1, sample_module)
							a2 = nu.sample_noise(a2, sample_module)

						corr = pearsonr(a1, a2)
						#print corr
						corr_sum += corr[0]
						#print corr_sum


					D[i,j] = corr_sum/average_cnt
					#print "D %s" % (D[i,j],)

				else:
					corr = pearsonr(np.array(data_set[i,])[0], np.array(data_set[j,])[0])
					D[i,j] = corr[0]

			elif dependence_metric == "mutual_information":
				a_i = np.array(data_set[i,])[0]
				a_j = np.array(data_set[j,])[0]
				a_i_list = [[x] for x in a_i.tolist()]
				a_j_list = [[x] for x in a_j.tolist()]
				minfo= ee.mi(a_i_list, a_j_list)
				D[i,j] = minfo
			elif dependence_metric == "mutual_information_score":
				a_i = np.array(data_set[i,])[0]
				a_j = np.array(data_set[j,])[0]
				minfo= skm.mutual_info_score(a_i, a_j)
				D[i,j] = minfo

	Ddist = []
	for i, j in it.combinations(range(data_set.shape[0]),2):
		if dependence_metric == "pearson":
			Ddist.append( 1 - D[i,j] )
		elif dependence_metric == "mutual_information" or dependence_metric == "mutual_information_score":
			#kdrew: this makes an actual distance, not sure of the right way of doing this
			Ddist.append( D.max() - D[i,j] )
			#kdrew: appeared random
			#Ddist.append( D[i,i] - D[i,j] )


	D = np.nan_to_num(D)
	Ddist = np.nan_to_num(Ddist)
	#print D
	#Y = sch.linkage(D, method='centroid')
	Y = sch.linkage(Ddist, method=cluster_method)
	#Y2 = sch.linkage(D, method='complete')
	#Y2 = sch.linkage(D, method='single')
	Y2 = sch.linkage(Ddist, method=cluster_method)

	return Y, Y2, D

def plotJustDendrogram(Y, plot_filename, new_id_map, circle_annotate=None, title=""):
	fig = pylab.figure(figsize=(14,14))
	#ax1 = fig.add_axes([0.09, 0.1, 0.11, 0.6])
	ax1 = pylab.subplot2grid((len(new_id_map),14),(0,0), rowspan=len(new_id_map), colspan=14)
	#pylab.title(title)
	dendrogram = sch.dendrogram(Y, orientation='right')
	#dendrogram = sch.dendrogram( Y )
	print dendrogram['leaves']
	print [new_id_map[z] for z in dendrogram['leaves']]

	#kdrew: testing purposes
	#circle_annotate = range(len(dendrogram['icoord']))
	if circle_annotate != None:
		for j, k in enumerate(zip(dendrogram['icoord'], dendrogram['dcoord'])):
			i = k[0]
			d = k[1]
			print j
			print d[1]
			cl_id = Y[:,2].tolist().index(d[1])
			x = 0.5 * sum(i[1:3])
			y = d[1]
			ax1.scatter(y, x, s=circle_annotate[cl_id], alpha=0.5)
			#ax1.annotate(str(circle_annotate[cl_id]), (y,x))
			#ax1.annotate(str(cl_id), (y,x))


	ax1.set_yticklabels([new_id_map[z] for z in dendrogram['leaves']])

	fig.savefig(plot_filename)
	pylab.close(fig)

def plotDendrogramProfile(data_set, Y, plot_filename, new_id_map, total_occupancy=False, title=""):

	if total_occupancy:
		col_sum = np.sum(data_set,0)
		shared_cols =  np.array(np.nonzero(col_sum)[1]).reshape(-1).tolist()
		print shared_cols
		data_set = data_set[np.ix_(range(0,data_set.shape[0]), shared_cols)]

	fig = pylab.figure(figsize=(14,14))
	#ax1 = fig.add_axes([0.09, 0.1, 0.11, 0.6])
	ax1 = pylab.subplot2grid((len(new_id_map),14),(0,0), rowspan=len(new_id_map), colspan=3)
	pylab.title(title)
	dendrogram = sch.dendrogram(Y, orientation='right')
	print dendrogram['leaves']
	print [new_id_map[z] for z in dendrogram['leaves']]
	ax1.set_yticklabels([new_id_map[z] for z in dendrogram['leaves']])

	data_subplots = []
	for i in xrange(len(dendrogram['leaves'])):
		ax = pylab.subplot2grid((len(new_id_map),14),(i,4), colspan=10)
		data_subplots.append(ax)

	max_value = np.max(data_set)
		
	for i, leaf in enumerate(dendrogram['leaves']):
		#kdrew: bottom of dendrogram is index 0, top is len(dendrogram['leaves'])
		j = len(dendrogram['leaves']) - 1 - i
		barcolor = "blue"
		#print "i: %s leaf: %s" % (i, leaf)
		data_row = data_set[leaf]
		#print data_row
		data_array_cols = np.array(data_row.reshape(-1))[0]
		#print data_array_cols
		#print data_array_cols[data_array_cols != 0]

		data_subplots[j].bar(np.arange(len(data_array_cols)), map(float,data_array_cols), align='center', facecolor=barcolor, alpha=0.5 )
		data_subplots[j].set_xlim(0, len(data_array_cols))
		data_subplots[j].axes.set_yticks([])
		#kdrew: only show tick marks on bottom profile 
		if i != 0:
			data_subplots[j].axes.set_xticks([])

	fig.savefig(plot_filename)
	pylab.close(fig)




#kdrew: code modified from stackoverflow
#http://stackoverflow.com/a/3011894
def plotDendrogram(Y, Y2, D, plot_filename, new_id_map, title=""):

	fig = pylab.figure(figsize=(8,8))
	ax1 = fig.add_axes([0.09, 0.1, 0.11, 0.6])
	dendrogram = sch.dendrogram(Y, orientation='right')
	print dendrogram['leaves']
	print [new_id_map[z] for z in dendrogram['leaves']]
	ax1.set_yticklabels([new_id_map[z] for z in dendrogram['leaves']], size='x-small')

	print new_id_map

	ax2 = fig.add_axes([0.3, 0.78, 0.6, 0.2])
	dendrogram2 = sch.dendrogram(Y2)
	ax2.set_xticklabels([new_id_map[z] for z in dendrogram2['leaves']], rotation='vertical', size='x-small')
	#pylab.title(title)

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
	pylab.close(fig)

	#print model.rows_
	#print model.columns_
	#print new_id_map

	return D


if __name__ == "__main__":
	main()


