

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
	parser.add_argument("--input_benchmark_pickle", action="store", dest="benchmark_filename", required=True, 
						help="Filename of pdb benchmark pickle")
	parser.add_argument("--data_threshold", action="store", dest="data_threshold", type=int, required=False, default=10,
						help="Only consider scores with data points above threshold")

	args = parser.parse_args()


	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.build_pdb_benchmark.py
	#mscpdbs = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/ms_complete_pdbs.p", "rb"))
	mscpdbs = cPickle.load(open(args.benchmark_filename, "rb"))

	#kdrew: pickle comes from protein_complex_maps.stoichiometry.benchmark.compute_benchmark_probs.py
	#mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_ms_complete_pdbs_results.p", "rb"))
	#mscpdbs_results = cPickle.load(open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.p", "rb"))
	mscpdbs_results = cPickle.load(open(args.results_filename, "rb"))

	#print mscpdbs_results
	#print mscpdbs

	score = 0
	top_ten = 0
	top_one = 0
	count = 0
	for pdb_id in mscpdbs_results:
		sorted_results = sorted(mscpdbs_results[pdb_id][0].iteritems(), key=operator.itemgetter(1))
		sorted_results.reverse()
		if mscpdbs_results[pdb_id][1] > args.data_threshold:
			count += 1
			print pdb_id
			print "true: %s" % (mscpdbs[pdb_id],)
			print "lengths: %s" % (pu.get_length_uniprot( mscpdbs[pdb_id].keys() ))
			print "data_points: %s" % (mscpdbs_results[pdb_id][1],)
			print "mapping: %s" % (mscpdbs_results[pdb_id][2],)
			true_map = {}
			for chain in mscpdbs_results[pdb_id][2]:
				true_map[chain] = mscpdbs[pdb_id][mscpdbs_results[pdb_id][2][chain]]
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
			
	print "score: %s" % (score,)
	print "top_ten: %s" % (top_ten,)
	print "top_one: %s" % (top_one,)
	print "count: %s" % (count,)



	#i = 0
	#for pdb_id in mscpdbs_results:
	#	if mscpdbs_results[pdb_id][1] > args.data_threshold:
	#		print "id: %s, ND: %s, meanRatio: %s, %s" % (pdb_id, top_rank_listND[i], top_rank_list[i], top_rank_listND[i] > top_rank_list[i])
	#		i+=1


if __name__ == "__main__":
	main()


