
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import itertools as it
import cPickle 
import operator
import ast

import protein_complex_maps.protein_util as pu


#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.build_pdb_benchmark.py
mscpdbs = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/ms_complete_pdbs.p", "rb"))

#kdrew: pickle comes from protein_complex_maps.stoichiometry.benchmark.compute_benchmark_probs.py
#mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_ms_complete_pdbs_results.p", "rb"))

mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.p", "rb"))

mscpdbs_noData_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_lenNormal_noData_ms_complete_pdbs_results.p", "rb"))

#plot_filename = "./HS_ms2_elutions_msds_lenNormal_noData_ms_complete_pdbs_results.pdf"
plot_filename = "./HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.pdf"
fig = plt.figure()
ax1 = fig.add_subplot(111)


#print mscpdbs_results
#print mscpdbs

def analyze( results ):
	top_rank_list = []
	score = 0
	top_ten = 0
	top_one = 0
	count = 0
	for pdb_id in results:
		sorted_results = sorted(results[pdb_id][0].iteritems(), key=operator.itemgetter(1))
		sorted_results.reverse()
		if results[pdb_id][1] > 10:
			count += 1
			print pdb_id
			print "true: %s" % (mscpdbs[pdb_id],)
			print "lengths: %s" % (pu.get_length_uniprot( mscpdbs[pdb_id].keys() ))
			print "data_points: %s" % (results[pdb_id][1],)
			print "mapping: %s" % (results[pdb_id][2],)
			true_map = {}
			for chain in results[pdb_id][2]:
				true_map[chain] = mscpdbs[pdb_id][results[pdb_id][2][chain]]
			print "true_map: %s" % (true_map,)
			top_rank = None
			for rank, pred in enumerate(sorted_results):
				stoich = ast.literal_eval(pred[0])
				if rank < 10:
					print "rank: %s, stoich: %s, score: %s" % (rank, stoich, pred[1])
				matched_all = True
				for comb in it.combinations(stoich.keys(),2):
					true_ratio = 1.0*true_map[comb[0]]/true_map[comb[1]]
					pred_ratio = 1.0*stoich[comb[0]]/stoich[comb[1]]
					if true_ratio != pred_ratio:
						matched_all = False

				#print "match?: %s" % (matched_all,)
				if matched_all and top_rank == None:
					top_rank  = rank

			print "top_rank: %s" % (top_rank,)
			score += top_rank
			if top_rank < 10:
				top_ten += 1
			if top_rank == 0:
				top_one += 1

			top_rank_list.append(top_rank)

	return (score, top_ten, top_one, count, top_rank_list)

score, top_ten, top_one, count, top_rank_list = analyze( mscpdbs_results )
scoreND, top_tenND, top_oneND, countND, top_rank_listND = analyze( mscpdbs_noData_results )
		
print "score: %s" % (score,)
print "top_ten: %s" % (top_ten,)
print "top_one: %s" % (top_one,)
print "count: %s" % (count,)

print "scoreND: %s" % (scoreND,)
print "top_tenND: %s" % (top_tenND,)
print "top_oneND: %s" % (top_oneND,)
print "countND: %s" % (countND,)


ind = np.arange(len(top_rank_list))
width = 0.35
#kdrew: this plots the top rank (prediction) of each pdb in benchmark
ax1.bar(ind+width, top_rank_listND, color='b', alpha=0.5)
ax1.bar(ind+width, top_rank_list, color='g', alpha=0.5)
ax1.set_xticks(ind+width)
ax1.set_xticklabels( ind, rotation=90, fontsize=6 )
plt.title("HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results")

plt.savefig(plot_filename)

