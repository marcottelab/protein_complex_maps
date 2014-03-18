
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
from random import shuffle

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
	parser.add_argument("--randomize", action="store_true", dest="randomize", required=False, default=False, 
						help="Flag for randomizing results of second pickle")
	parser.add_argument("--trials", action="store", type=int, dest="trials", required=False, default=1, 
						help="Number of randomized trials")
	parser.add_argument("--sort_mean", action="store_true", dest="sort_mean", required=False, default=False, 
						help="Flag to sort benchmark on mean")
	parser.add_argument("--sort_rank", action="store_true", dest="sort_rank", required=False, default=False, 
						help="Flag to sort benchmark on results1 rank")

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

	score, top_ten, top_one, count, top_rank_list, results_list, pdb_list = analyze( mscpdbs_results, mscpdbs, args )
	for result in results_list:
		print "%s\t%s" % (result[0], result[1])

 	top_rank_list = zip(*top_rank_list)[0]

	ind = np.arange(len(top_rank_list))
	width = 0.5

	if args.results_filename2 != None:
		mscpdbs_results2 = cPickle.load(open(args.results_filename2, "rb"))
		score2, top_ten2, top_one2, count2, top_rank_list2, results_list2, pdb_list2 = analyze( mscpdbs_results2, mscpdbs, args, randomize=args.randomize, trials=args.trials )

		assert(len(pdb_list) == len(pdb_list2))

		results.append(mscpdbs_results2)
		if args.trials > 1:     
			trial_means = []
			trial_std = []
			for i in top_rank_list2:
				trial_means.append(np.array(i).mean())
				trial_std.append(np.array(i).std())
		
			if args.sort_mean:
				trial_means, trial_std, top_rank_list, pdb_list = zip(*sorted(zip(trial_means, trial_std, top_rank_list, pdb_list)))
			if args.sort_rank:
				top_rank_list, trial_means, trial_std, pdb_list = zip(*sorted(zip(top_rank_list, trial_means, trial_std, pdb_list)))


			
			#ax1.errorbar(ind+width, trial_means, yerr=trial_std, fmt='o', color='b', alpha=0.5, ms=2)
			ax1.plot(ind+width, trial_means, '-o', color='b', alpha=0.5, ms=2)
			#ax1.plot(ind+width, np.array(trial_means)+np.array(trial_std), '--', color='b', alpha=0.5, ms=2)
			#ax1.plot(ind+width, np.array(trial_means)-np.array(trial_std), '--', color='b', alpha=0.5, ms=2)
			ax1.fill_between(ind+width, np.array(trial_means)+np.array(trial_std), np.array(trial_means)-np.array(trial_std), color='b', alpha=0.1)
		else:
			top_rank_list2 = zip(*top_rank_list2)[0]
			ax1.plot(ind+width, top_rank_list2, 'o', color='b', alpha=0.5, ms=2)

			for result in results_list2:
				print "%s\t%s" % (result[0], result[1])


	if args.plot_filename != None:


		#kdrew: this plots the top rank (prediction) of each pdb in benchmark
		ax1.plot(ind+width, top_rank_list, '-o', color='g', alpha=0.5, ms=2)
		ax1.set_xticks(ind+width)
		stoich_list = [mscpdbs[x].values() for x in pdb_list]
		#ax1.set_xticklabels( ind, rotation=90, fontsize=6 )
		ax1.set_xticklabels( stoich_list, rotation=90, fontsize=4 )
		#plt.title("HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results")

		ylim = plt.ylim()
		plt.ylim((-3, ylim[1]))

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
		for pdb_id in pdb_list:
			if min(mscpdbs_results[pdb_id][1]) > args.data_threshold:
				#print "id: %s, Rank2: %s, Rank1: %s, %s" % (pdb_id, np.array(top_rank_list2[i]).mean(), top_rank_list[i], np.array(top_rank_list2[i]).mean() > top_rank_list[i])
				print "id: %s, Rank2: %s, Rank1: %s, TrueStoich: %s " % (pdb_id, trial_means[i], top_rank_list[i], mscpdbs[pdb_id])
				i+=1


		z = (np.array(score2).mean() - score)/np.array(score2).std()
		print "z(score): %s" % (z,)

#print mscpdbs_results
#print mscpdbs

def analyze( results, mscpdbs, args, randomize=False, trials=1):
	results_list = []
	top_rank_list = []
	pdb_list = []

	score_list = []
	top_ten_list = []
	top_one_list = []
	count_list = []

	for i in xrange(trials):

		score = 0
		top_ten = 0
		top_one = 0
		count = 0

		tmp_top_rank_list = []
		for pdb_id in results:
			sorted_results = sorted(results[pdb_id][0].iteritems(), key=operator.itemgetter(1))
			if randomize:
				#kdrew: shuffle results and resort
				shuffle(sorted_results)
			else:
				sorted_results = sorted(sorted_results, key=lambda tup: tup[1])
				sorted_results.reverse()

			if min(results[pdb_id][1]) > args.data_threshold:

				count += 1
				print pdb_id
				print "true: %s" % (mscpdbs[pdb_id],)
				#print "lengths: %s" % (pu.get_length_uniprot( mscpdbs[pdb_id].keys() ))

				print "data_points: %s" % (results[pdb_id][1],)
				print "mapping: %s" % (results[pdb_id][2],)
				true_map = {}
				for chain in results[pdb_id][2]:
					true_map[chain] = mscpdbs[pdb_id][results[pdb_id][2][chain]]
				
				print "true_map: %s" % (true_map,)
				top_rank = None
				for rank, pred in enumerate(sorted_results):
					stoich = ast.literal_eval(pred[0])

					#if rank < 10:
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
						results_list.append((1, rank))
					else:
						results_list.append((-1, rank))

				print "top_rank: %s" % (top_rank,)
				score += top_rank
				if top_rank < 10:
					top_ten += 1
				if top_rank == 0:
					top_one += 1

				tmp_top_rank_list.append(top_rank)

				if i < 1:
					pdb_list.append(pdb_id)

		if len(top_rank_list) == 0:
			#top_rank_list = tmp_top_rank_list
			top_rank_list = map(list,zip(tmp_top_rank_list))
		else:
			for j in xrange(len(top_rank_list)):
				top_rank_list[j].append(tmp_top_rank_list[j])

		score_list.append(score)
		top_ten_list.append(top_ten)
		top_one_list.append(top_one)
		count_list.append(count)

	#if trials > 1:
	#	ret_top_rank_list = [1.0*x/trials for x in top_rank_list]
	#	score = 1.0*score/trials
	#	top_ten = 1.0*top_ten/trials
	#	top_one = 1.0*top_one/trials
	#	count = 1.0*count/trials
	#else:

	return (score_list, top_ten_list, top_one_list, count_list, top_rank_list, results_list, pdb_list)



if __name__ == "__main__":
	main()

