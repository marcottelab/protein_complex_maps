
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import itertools as it
import cPickle 
import operator
import ast
import argparse

import protein_complex_maps.protein_util as pu

def main():
	parser = argparse.ArgumentParser(description="Tool to analyze pdb stoichometry benchmark results")
	parser.add_argument("--results_pickle", action="store", dest="results_filename", required=True, 
						help="Filename of compute_benchmark_probs pickle")
	parser.add_argument("--results_pickle2", action="store", dest="results_filename2", required=False, default=None, 
						help="Filename of compute_benchmark_probs pickle for comparison")
	parser.add_argument("--input_benchmark_pickle", action="store", dest="benchmark_filename", required=True, 
						help="Filename of pdb benchmark pickle")
	parser.add_argument("--data_threshold", action="store", dest="data_threshold", type=int, required=False, default=10,
						help="Only consider scores with data points above threshold")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None, 
						help="Filename of results plot")

	args = parser.parse_args()


	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.build_pdb_benchmark.py
	#mscpdbs = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/ms_complete_pdbs.p", "rb"))
	mscpdbs = cPickle.load(open(args.benchmark_filename, "rb"))

	#kdrew: pickle comes from protein_complex_maps.stoichiometry.benchmark.compute_benchmark_probs.py
	#mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_ms_complete_pdbs_results.p", "rb"))
	#mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.p", "rb"))
	mscpdbs_results = cPickle.load(open(args.results_filename, "rb"))

	results = [mscpdbs_results]

	#print mscpdbs_results
	#print mscpdbs

	#plot_filename = "./HS_ms2_elutions_msds_lenNormal_noData_ms_complete_pdbs_results.pdf"
	#plot_filename = "./HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.pdf"
	fig = plt.figure()
	ax1 = fig.add_subplot(111)

	score, top_ten, top_one, count, top_rank_list = analyze( mscpdbs_results, mscpdbs, args )

	ind = np.arange(len(top_rank_list))
	width = 0.35

	if args.results_filename2 != None:
		mscpdbs_results2 = cPickle.load(open(args.results_filename2, "rb"))
		score2, top_ten2, top_one2, count2, top_rank_list2 = analyze( mscpdbs_results2, mscpdbs, args )
		results.append(mscpdbs_results2)
		ax1.bar(ind+width, top_rank_list2, color='b', alpha=0.5)

	if args.plot_filename != None:

		#kdrew: this plots the top rank (prediction) of each pdb in benchmark
		ax1.bar(ind+width, top_rank_list, color='g', alpha=0.5)
		ax1.set_xticks(ind+width)
		ax1.set_xticklabels( ind, rotation=90, fontsize=6 )
		#plt.title("HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results")

		plt.savefig(args.plot_filename)

			
	print "score: %s" % (score,)
	print "top_ten: %s" % (top_ten,)
	print "top_one: %s" % (top_one,)
	print "count: %s" % (count,)

	if args.results_filename2 != None:
		print "score2: %s" % (score2,)
		print "top_ten2: %s" % (top_ten2,)
		print "top_one2: %s" % (top_one2,)
		print "count2: %s" % (count2,)


		i = 0
		for pdb_id in mscpdbs_results:
			if mscpdbs_results[pdb_id][1] > args.data_threshold:
				print "id: %s, 2: %s, meanRatio: %s, %s" % (pdb_id, top_rank_list2[i], top_rank_list[i], top_rank_list2[i] > top_rank_list[i])
				i+=1

#print mscpdbs_results
#print mscpdbs

def analyze( results, mscpdbs, args ):
	top_rank_list = []

	score = 0
	top_ten = 0
	top_one = 0
	count = 0

	for pdb_id in results:
		sorted_results = sorted(results[pdb_id][0].iteritems(), key=operator.itemgetter(1))
		sorted_results.reverse()
		if results[pdb_id][1] > args.data_threshold:

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



if __name__ == "__main__":
	main()

