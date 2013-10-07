
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



def plot_bicluster(data_set, bicluster1):
	data_subplots = []
	f, data_subplots = plt.subplots(len(data_set),1,sharex='col')

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

		data_array_in_cols[np.ix_(list(cols - set(bicluster1.columns())))] = 0

		data_subplots[i].bar(np.arange(len(data_array_in_cols)), map(float,data_array_in_cols), align='center', facecolor=in_barcolor, alpha=0.5 )
		data_subplots[i].bar(np.arange(len(data_array_out_cols)), map(float,data_array_out_cols), align='center', facecolor=out_barcolor, alpha=0.5 )

		data_subplots[i].axes.set_yticklabels([],visible=False)
		data_subplots[i].set_ylabel(i,rotation='horizontal', color=in_barcolor)

	plt.show()
