

import glob
import cPickle
import pickle

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


def main():

	msds = None
	length_normalize = False
	threshold_normalize = False

	#kdrew: this will read in the raw msds pickle and normalize by length as well as map ids
	if length_normalize:
		#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
		msds = pickle.load( open( "./HS_ms2_elutions_msds.p", "rb" ) )
		print msds.get_id_dict()
		#kdrew: map ids
		msds.map_ids("ENSEMBL_ID", "ACC")
		print msds.get_id_dict()

		#kdrew: normalize length
		dm_length = nu.normalize_length( msds.get_data_matrix(), msds.get_id_dict() )
		msds.set_data_matrix(dm_length)

		print msds.get_data_matrix()

		cPickle.dump( msds, open("./HS_ms2_elutions_msds_lenNormal.p","wb"))

	#kdrew: length normalization is already done, just read in pickle
	else:
		#kdrew: pickle comes from running protein_complex_maps.stoichiometry.benchmark.compute_benchmark_probs.py with length_normalize = True
		msds = pickle.load( open( "./HS_ms2_elutions_msds_lenNormal.p", "rb" ) )
		print msds.get_id_dict()
	
	if threshold_normalize:
		msds_dm = msds.get_data_matrix()
		msds_dm = nu.threshold(msds_dm, threshold=0.0001)
		msds.set_data_matrix(msds_dm)

		cPickle.dump( msds, open("./HS_ms2_elutions_msds_lenNormal_threshold.0001.p","wb"))


if __name__ == "__main__":
	main()

