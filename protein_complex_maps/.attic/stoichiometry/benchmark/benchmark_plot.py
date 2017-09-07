
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cPickle 
import operator


plot_filename = "./ms_complete_pdbs_benchmark.pdf"

mscpdbs = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/ms_complete_pdbs.p", "rb"))

mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.p", "rb"))

fig = plt.figure()
ax1 = fig.add_subplot(111)

stoich_count = dict()
have_data_count = dict()
for pdb_id in mscpdbs:
	print mscpdbs[pdb_id]
	stoich = mscpdbs[pdb_id].values()
	stoich.sort()
	try:
		stoich_count[tuple(stoich)] +=1
	except KeyError:
		stoich_count[tuple(stoich)] = 1

for pdb_id in mscpdbs:
	stoich = mscpdbs[pdb_id].values()
	stoich.sort()
	try:
		if pdb_id in mscpdbs_results.keys() and mscpdbs_results[pdb_id][1] > 10:
			have_data_count[tuple(stoich)] +=1
	except KeyError:
		if pdb_id in mscpdbs_results.keys() and mscpdbs_results[pdb_id][1] > 10:
			have_data_count[tuple(stoich)] = 1

stoich_count_list = sorted(stoich_count.iteritems(), key=operator.itemgetter(1))
stoich_count_list.reverse()
print "stoich_count_list: %s" % (stoich_count_list,)

print "have_data_count: %s" % (have_data_count,)

have_data_count_list = []
for stoich in zip(*stoich_count_list)[0]:
	try:
		have_data_count_list.append(have_data_count[stoich])
	except KeyError:
		have_data_count_list.append(0)

ind = np.arange(len(stoich_count_list))
width = 0.35
#kdrew: this plots the number of pdbs in full benchmark by stoichiometry
ax1.bar(ind+width, zip(*stoich_count_list)[1])
#kdrew: this plots the number of pdbs in full benchmark by stoichiometry that we have > 10 data points for
ax1.bar(ind+width, have_data_count_list, color='g')
ax1.set_xticks(ind+width)
ax1.set_xticklabels( zip(*stoich_count_list)[0], rotation=45, fontsize=8 )
plt.title("PDB Stoichiometry Benchmark")

plt.savefig(plot_filename)

