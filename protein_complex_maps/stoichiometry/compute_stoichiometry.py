import glob
import cPickle
import pickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


def main():

	parser = argparse.ArgumentParser(description="Tool to compute stoichiometry probabilities (likelihoods)")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle")
	parser.add_argument("--input_complex_list", action="store", dest="complex_filename", required=True, 
						help="Filename of complex protein identifiers")
	parser.add_argument("--stoichiometry_filename", action="store", dest="stoichiometry_filename", required=True, 
						help="Filename of stoichiometries from pdb website (slightly massaged)")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
						help="Filename of output pickle")
	parser.add_argument("--prior_type", action="store", dest="prior_type", required=False, default="uniform",
						help="Type of prior to use for calculation, [uniform, pdb]")

	args = parser.parse_args()

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

	complex_file = open(args.complex_filename, 'rb')
	complex_list = []
	for i in complex_file.readlines():
		complex_list.append(i.split())
	


	for prot_ids in complex_list:
		if msds_id_set >= set(prot_ids):
			print "complex ids are in msds"
			results, num_dpts, prot_ids_map = rs.relative_stoichiometry( msds, prot_ids, s_set, prior_type=args.prior_type )
			print results
			results_dict[prot_ids.__str__()] = (results, num_dpts, prot_ids_map)
		else:
			print "complex ids are NOT in msds"

	cPickle.dump(results_dict, open(args.output_filename, "wb"))

if __name__ == "__main__":
	main()

