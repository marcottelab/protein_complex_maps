

import glob
import cPickle
import pickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


def main():

	parser = argparse.ArgumentParser(description="Tool to compute stoichiometry probabilities (likelihoods) for a pdb benchmark")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle")
	parser.add_argument("--input_benchmark_pickle", action="store", dest="benchmark_filename", required=True, 
						help="Filename of pdb benchmark pickle")
	parser.add_argument("--pdb_ids", action="store", dest="pdb_ids", required=False, nargs='+', default=None,
						help="List of pdbs in benchmark to compute")
	parser.add_argument("--stoichiometry_filename", action="store", dest="stoichiometry_filename", required=True, 
						help="Filename of stoichiometries from pdb website (slightly massaged)")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
						help="Filename of output pickle")
	parser.add_argument("--prior_type", action="store", dest="prior_type", required=False, default="uniform",
						help="Type of prior to use for calculation, [uniform, pdb]")

	args = parser.parse_args()

	#kdrew: read in benchmark
	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.build_pdb_benchmark.py
	ms_complete_pdbs = cPickle.load(open(args.benchmark_filename, "rb"))

	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.normalize_benchmark_msds.py with length_normalize = True
	msds = pickle.load( open( args.msds_filename, "rb" ) )
	print msds.get_id_dict()


	#kdrew: read in stoichiometries, file is from pdb website (slightly massaged)
	#filename = "/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/pdb_stoichiometry_list.txt"
	sample_file = open(args.stoichiometry_filename, 'rb')
	s_set = st.Stoichiometries()
	s_set.read_stoichiometries(sample_file)

	msds_id_dict = msds.get_id_dict()
	msds_id_set = set(msds_id_dict.keys())
	results_dict = dict()

	if args.pdb_ids == None:
		args.pdb_ids = ms_complete_pdbs.keys()
	#kdrew: for every benchmark complex, call relative stoichiometry
	for pdb_id in args.pdb_ids:
		print pdb_id
		prot_ids = ms_complete_pdbs[pdb_id].keys()
		if msds_id_set >= set(prot_ids):
			print "benchmark ids are in msds"
			results, num_dpts, prot_ids_map, single_class_flag = rs.relative_stoichiometry( msds, prot_ids, s_set, prior_type=args.prior_type )
			print "results: %s" % (results,)
			if single_class_flag:
				results_dict[pdb_id] = (results, num_dpts, prot_ids_map)
			else:
				print "GMM did not converge to single class on all cases, single_class_flag == False"
		else:
			print "benchmark ids are NOT in msds"

	cPickle.dump(results_dict, open(args.output_filename, "wb"))

if __name__ == "__main__":
	main()

