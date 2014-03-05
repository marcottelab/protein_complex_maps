
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import multiprocessing as mp
import argparse
import cPickle
import pickle

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

	args = parser.parse_args()

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )

	plot_profile(msds, args.proteins, total_occupancy=args.total_occupancy, savefilename=args.plot_filename)

def plot_profile(msds, protein_ids, total_occupancy=False, ylim_max=False, savefilename=None):
	data_subplots = []
	data_set, new_id_map = msds.get_subdata_matrix(protein_ids) 

	#kdrew: only look at columns that have all genes with values
	if total_occupancy:
		col_sum = np.sum(data_set,0)
		shared_cols =  np.array(np.nonzero(col_sum)[1]).reshape(-1).tolist()
		print shared_cols
		data_set = data_set[np.ix_(range(0,data_set.shape[0]), shared_cols)]

	print data_set

	f, data_subplots = plt.subplots(len(data_set),1,sharex='col')

	max_value = np.max(data_set)

	for i, data_row in enumerate(data_set):
		barcolor = "blue"

		data_array_cols = np.array(data_row.reshape(-1))[0]

		cols = set(range(0,len(data_array_cols)))

		data_subplots[i].bar(np.arange(len(data_array_cols)), map(float,data_array_cols), align='center', facecolor=barcolor, alpha=0.5 )

		#data_subplots[i].axes.set_yticklabels([],visible=False)

		if ylim_max:
			data_subplots[i].axes.set_ylim(0,max_value)

		data_subplots[i].set_ylabel(new_id_map[i],rotation='horizontal', color=barcolor)

	if savefilename is None:
		plt.show()
	else:
		plt.savefig(savefilename)


if __name__ == "__main__":
	main()



