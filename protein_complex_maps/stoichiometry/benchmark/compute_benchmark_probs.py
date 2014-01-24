

import glob
import cPickle
import pickle

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


def main():

	#kdrew: read in benchmark
	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.build_pdb_benchmark.py
	ms_complete_pdbs = cPickle.load(open("./ms_complete_pdbs.p", "rb"))

	msds = None
	#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.normalize_benchmark_msds.py with length_normalize = True
	msds = pickle.load( open( "./HS_ms2_elutions_msds_lenNormal.p", "rb" ) )
	print msds.get_id_dict()


	#kdrew: read in stoichiometries, file is from pdb website (slightly massaged)
	filename = "/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/pdb_stoichiometry_list.txt"
	sample_file = open(filename, 'rb')
	s_set = st.Stoichiometries()
	s_set.read_stoichiometries(sample_file)

	msds_id_dict = msds.get_id_dict()
	msds_id_set = set(msds_id_dict.keys())
	results_dict = dict()
	#kdrew: for every benchmark complex, call relative stoichiometry
	for pdb_id in ms_complete_pdbs:
		print pdb_id
		prot_ids = ms_complete_pdbs[pdb_id].keys()
		if msds_id_set >= set(prot_ids):
			print "benchmark ids are in msds"
			results, num_dpts, prot_ids_map = rs.relative_stoichiometry( msds, prot_ids, s_set )
			print results
			results_dict[pdb_id] = (results, num_dpts, prot_ids_map)
		else:
			print "benchmark ids are NOT in msds"

	cPickle.dump(results_dict, open("./HS_ms2_elutions_msds_lenNormal_mean_ratio_true_ms_complete_pdbs_results.p", "wb"))

if __name__ == "__main__":
	main()

