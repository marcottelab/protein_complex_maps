
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import multiprocessing as mp


def plot_score_distribution(distribution_dict, score=None, savefilename=None):
	p = mp.Process(target=plot_score_distribution_worker, args=(distribution_dict, score, savefilename))
	p.start()

def plot_score_distribution_worker(distribution_dict, score=None, savefilename=None):
	f, data_subplots = plt.subplots(len(distribution_dict.keys()),1)
	for i, t in enumerate(distribution_dict.keys()):
		data_subplots[i].hist( distribution_dict[t], bins=50 )
		data_subplots[i].set_ylabel(t,rotation='horizontal')

	if savefilename is None:
		plt.show()
	else:
		plt.savefig(savefilename)


#kdrew: wrapper function which spawns processor so plotting does not slow up driver script
def plot_bicluster(data_set, bicluster1, ylim_max=False, savefilename=None):
	p = mp.Process(target=plot_bicluster_worker, args=(data_set, bicluster1, ylim_max, savefilename))
	p.start()

def plot_bicluster_worker(data_set, bicluster1, ylim_max=False, savefilename=None):
	data_subplots = []
	f, data_subplots = plt.subplots(len(data_set),1,sharex='col')

	max_value = np.max(data_set)

	for i, data_row in enumerate(data_set):
		in_barcolor = "0.5"
		out_barcolor = "1.0"
		if i in bicluster1.rows():
			in_barcolor = "blue"
			out_barcolor = "1.0"

		data_array_in_cols = np.array(data_row.reshape(-1))[0]
		data_array_out_cols = np.array(data_row.reshape(-1))[0]

		cols = set(range(0,len(data_array_in_cols)))

		data_array_out_cols[np.ix_(bicluster1.columns())] = 0

		if len(list(cols - set(bicluster1.columns()))) != 0:
			data_array_in_cols[np.ix_(list(cols - set(bicluster1.columns())))] = 0

		data_subplots[i].bar(np.arange(len(data_array_in_cols)), map(float,data_array_in_cols), align='center', facecolor=in_barcolor, alpha=0.5 )
		data_subplots[i].bar(np.arange(len(data_array_out_cols)), map(float,data_array_out_cols), align='center', facecolor=out_barcolor, alpha=0.5 )

		data_subplots[i].axes.set_yticklabels([],visible=False)
		data_subplots[i].set_ylabel(i,rotation='horizontal', color=in_barcolor)

		if ylim_max:
			data_subplots[i].axes.set_ylim(0,max_value)

	if savefilename is None:
		plt.show()
	else:
		plt.savefig(savefilename)


