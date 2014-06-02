
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import multiprocessing as mp
import argparse
import cPickle
import pickle

import protein_complex_maps.protein_util as pu

def main():
	parser = argparse.ArgumentParser(description="Plot abundance profile of proteins")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
						help="Protein ids in which to plot")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--total_occupancy", action="store_true", dest="total_occupancy", required=False, default=False,
						help="Flag to only plot columns where every gene has greater than zero counts")
	parser.add_argument("--fraction_range", action="store", dest="fraction_range", nargs='+', required=False, 
						help="Sets the range of fractions to plot")
	parser.add_argument("--genenames", action="store_true", dest="genenames", required=False, default=False,
						help="Set labels to genenames")

	args = parser.parse_args()

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )

	if args.fraction_range:
		frange = [int(x) for x in args.fraction_range]
	else:
		frange = None

	plot_profile(msds, args.proteins, total_occupancy=args.total_occupancy, savefilename=args.plot_filename, fraction_range=frange, genenames=args.genenames)

def plot_profile(msds, protein_ids, total_occupancy=False, ylim_max=False, savefilename=None, fraction_range = None, genenames=False):
	data_set, new_id_map = msds.get_subdata_matrix(protein_ids) 

	if fraction_range != None:
		print fraction_range
		#data_set = data_set[np.ix_(range(0,data_set.shape[0]), range(fraction_range[0], fraction_range[1]))]


	#kdrew: only look at columns that have all genes with values, CAUTION: does not do this as implemented
	#kdrew: just gets rid of zeros
	if total_occupancy:
		col_sum = np.sum(data_set,0)
		shared_cols =  np.array(np.nonzero(col_sum)[1]).reshape(-1).tolist()
		print shared_cols
		data_set = data_set[np.ix_(range(0,data_set.shape[0]), shared_cols)]

	#print data_set

	plot_profile_dataset(data_set, new_id_map, ylim_max, savefilename, fraction_range)

def plot_profile_dataset(data_set, id_map, ylim_max=False, savefilename=None, fraction_range=None, x_highlight=(0,0), y_highlight=(0,0), genenames=False):
	data_subplots = []

	if genenames:
		genename_map = pu.get_genenames_uniprot( id_map.values() )

	f, data_subplots = plt.subplots(len(data_set),1,sharex='col')

	max_value = np.max(data_set)

	for i, data_row in enumerate(data_set):
		barcolor = "blue"

		data_array_cols = np.array(data_row.reshape(-1))[0]

		#cols = set(range(0,len(data_array_cols)))

		data_subplots[i].bar(np.arange(len(data_array_cols)), map(float,data_array_cols), align='center', facecolor=barcolor, alpha=0.5 )

		data_subplots[i].axvspan(y_highlight[0], y_highlight[1], color='yellow', alpha=0.2)

		if i in xrange(x_highlight[0], x_highlight[1]):
			data_subplots[i].axvspan(0, len(data_array_cols), color='blue', alpha=0.2)

		#data_subplots[i].axes.set_yticklabels([],visible=False)

		if ylim_max:
			data_subplots[i].axes.set_ylim(0,max_value)

		if fraction_range:
			data_subplots[i].axes.set_xlim(fraction_range[0], fraction_range[1])

		try:
			if genenames:
				data_subplots[i].set_ylabel(genename_map[id_map[i]],rotation='horizontal', color=barcolor, fontsize=10)
			else:
				data_subplots[i].set_ylabel(id_map[i],rotation='horizontal', color=barcolor, fontsize=10)
		except KeyError:
			data_subplots[i].set_ylabel(id_map[i],rotation='horizontal', color=barcolor, fontsize=10)


		#data_subplots[i].axes.set_yticks(data_subplots[i].axes.get_yticks()[0::5])
		data_subplots[i].axes.set_yticks([])

	if savefilename is None:
		print "savefilename is None"
		plt.show()
	else:
		plt.savefig(savefilename)


if __name__ == "__main__":
	main()



