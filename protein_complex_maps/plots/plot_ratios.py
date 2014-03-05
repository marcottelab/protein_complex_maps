
import itertools
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os.path
from sklearn import mixture

from multiprocessing import Pool
import multiprocessing as mp
import argparse
import cPickle
import pickle

def main():
	parser = argparse.ArgumentParser(description="Plots ratios of proteomic profiles for given proteins ")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
						help="Protein ids in which to plot")
	parser.add_argument("--input_complex_list", action="store", dest="complex_filename", required=False, 
						help="Filename of complex protein identifiers")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--total_occupancy", action="store_true", dest="total_occupancy", required=False, default=False,
						help="Flag to only plot columns where every gene has greater than zero counts")
	parser.add_argument("--log_ratio", action="store_true", dest="log_ratio", required=False, default=False,
						help="Flag to log transform ratios")
	parser.add_argument("--histogram", action="store_true", dest="histogram", required=False, default=False,
						help="Plot ratios in a histogram")
	parser.add_argument("--data_points_threshold", action="store", dest="threshold", required=False, default=10,
						help="Pairs of proteins must have atleast specified number of data points, default=10")
	parser.add_argument("-j", action="store", dest="numOfProcs", required=False, default=1,
						help="Number of processors to use, default=1")

	args = parser.parse_args()

	pool = Pool(processes=int(args.numOfProcs))
	print "setup pool"

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )
	msds_id_dict = msds.get_id_dict()
	msds_id_set = set(msds_id_dict.keys())

	if args.complex_filename != None:

		complex_file = open(args.complex_filename, 'rb')
		complex_list = []
		for i in complex_file.readlines():
			complex_list.append(i.split())

		for i, prot_ids in enumerate(complex_list):
			if msds_id_set >= set(prot_ids):                                                                              
				plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(i)+os.path.splitext(args.plot_filename)[1]
				print "sending job"
				pool.apply_async(plot_ratios, (msds, prot_ids), dict(total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, savefilename=plot_fname, histogram=args.histogram, data_threshold=args.threshold))
				print "sent job"

	elif args.proteins != None:
		plot_ratios(msds, args.proteins, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, savefilename=args.plot_filename, histogram=args.histogram, data_threshold=args.threshold)
	else:
		print "use either --proteins or --complex_filename to specify protein ids"

	pool.close()
	pool.join()

def plot_ratios(msds, protein_ids, total_occupancy=False, log_ratio=False, ylim_max=False, savefilename=None, histogram=False, data_threshold=10):

	print total_occupancy
	print log_ratio
	print savefilename
	print histogram
	print data_threshold
	data_subplots = []
	data_set, new_id_map = msds.get_subdata_matrix(protein_ids) 

	print protein_ids
	#kdrew: only look at columns that have all genes with values
	if total_occupancy:
		col_sum = np.sum(data_set,0)
		shared_cols =  np.array(np.nonzero(col_sum)[1]).reshape(-1).tolist()
		print "shared_cols: %s" % (shared_cols,)
		data_set = data_set[np.ix_(range(0,data_set.shape[0]), shared_cols)]

	print "printing dataset"
	print "data_set: %s" % (data_set,)

	#kdrew: calculate all (protein) v all ratios of peak intensities
	ratio_dict = dict()
	#kdrew: for every protein
	for i, data_r in enumerate(data_set):
		data_row = np.array(data_r.reshape(-1))[0]

		#kdrew: for every other protein
		for j, data_r2 in enumerate(data_set):
			if i == j:
				continue

			ratio_dict[i,j] = []
			data_row2 = np.array(data_r2.reshape(-1))[0]
			#kdrew: for every fraction in data_matrix
			for k in range(len(data_row2)):
				#kdrew: make sure both proteins were seen at this fraction, 
				#kdrew: this might be a problem where floats are sometimes not exactly 0.0
				if data_row[k] != 0.0 and data_row2[k] != 0.0:
					print "%s : %s : %s" % (data_row[k], data_row2[k], data_row[k]/data_row2[k])
					if log_ratio:
						#kdrew: add ratio to list
						ratio_dict[i,j].append(np.log(data_row[k]/data_row2[k]))
					else:
						ratio_dict[i,j].append(data_row[k]/data_row2[k])

			
	print ratio_dict

	plot_flag = False
	f, data_subplots = plt.subplots(len(ratio_dict.keys()),1,sharex='col')

	for i, pair in enumerate(ratio_dict):
		if len(ratio_dict[pair]) < data_threshold or np.isnan(np.sum(ratio_dict[pair])):
			continue

		if histogram:
			#data_subplots[i].hist(ratio_dict[pair])

			print ratio_dict[pair]
			print new_id_map[pair[0]], new_id_map[pair[1]]
			color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm'])
			clf = mixture.DPGMM(n_components=5, cvtype='diag')
 			X = np.array([[x,] for x in ratio_dict[pair]])
			clf.fit(X)
			print clf.converged_
			if clf.converged_:
				print clf

			Y = clf.predict(X)
			print Y
			for j, (mean, color) in enumerate(zip(clf.means, color_iter)):
				class_data = np.array(ratio_dict[pair])[Y == j]
				if len(class_data) > 0:
					data_subplots[i].hist(class_data, color=color)
					data_subplots[i].set_title("%s : %s, converged? %s" % (new_id_map[pair[0]], new_id_map[pair[1]], clf.converged_))
					plot_flag = True



		else:
			barcolor = "blue"
			meancolor = "red"
			mediancolor = "orange"

			#kdrew: add two extra cols for mean and median
			cols = len(ratio_dict[pair])+2

			#kdrew: this gets a little hacky to put mean and median on the end of the bar plot
			ratio_array = np.array(ratio_dict[pair])
			print "mean %s, std %s, median %s" % ( ratio_array.mean(), ratio_array.std(), np.median(ratio_array))
			mean_list = np.repeat(0.0,cols-2).tolist()
			mean_list.append(ratio_array.mean())
			mean_list.append(0.0)
			median_list = np.repeat(0.0,cols-1).tolist()
			median_list.append(np.median(ratio_array))
			ratio_dict[pair].append(0.0)
			ratio_dict[pair].append(0.0)

			data_subplots[i].bar(np.arange(cols), map(float,ratio_dict[pair]), align='center', facecolor=barcolor, alpha=0.5 )
			data_subplots[i].bar(np.arange(cols), mean_list, align='center', facecolor=meancolor, alpha=0.5 )
			data_subplots[i].bar(np.arange(cols), median_list, align='center', facecolor=mediancolor, alpha=0.5 )

			if ylim_max:
				data_subplots[i].axes.set_ylim(0,max_value)

			data_subplots[i].set_ylabel(pair,rotation='horizontal', color=barcolor)



	if savefilename is None:
		plt.show()
	elif plot_flag:
		plt.savefig(savefilename)


if __name__ == "__main__":
	main()



