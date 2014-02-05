

import glob
import cPickle
import pickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


def main():

	parser = argparse.ArgumentParser(description="Tool to normalize Mass Spec Data Set (MSDS) pickle")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle")
	parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
						help="Output filename of normalized MSDS pickle")
	parser.add_argument("--length_normalize", action="store_true", dest="length_normalize", required=False, default=False,
						help="Normalize by length, currently requires ENSEMBL_ID")
	parser.add_argument("--peptide_normalize", action="store_true", dest="peptide_normalize", required=False, default=False,
						help="Normalize by peptide, requires peptide_digest_file flag")
	parser.add_argument("--peptide_digest_filename", action="store", dest="peptide_digest_filename", required=False, 
						help="Filename of possible digested peptides")
	parser.add_argument("--peptide_digest_format", action="store", dest="peptide_digest_format", required=False, default="raw",
						help="Format of peptide digest file, either raw or counts")
	parser.add_argument("--threshold_normalize", action="store_true", dest="threshold_normalize", required=False, default=False,
						help="Normalize by threshold")
	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, 
						help="Threshold for which to normalize by")
	parser.add_argument("--map_ids", action="store_true", dest="map_ids", required=False, default=False,
						help="Map ENSEMBL_IDs to Uniprot accs")

	args = parser.parse_args()

	#length_normalize = False
	#threshold_normalize = False

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )
	if args.map_ids:
		msds.map_ids("ENSEMBL_ID", "ACC")

	#kdrew: normalize by length as well as map ids
	if args.length_normalize:
		print msds.get_id_dict()
		#kdrew: map ids
		msds.map_ids("ENSEMBL_ID", "ACC")
		print msds.get_id_dict()

		#kdrew: normalize length
		dm_length = nu.normalize_length( msds.get_data_matrix(), msds.get_id_dict() )
		msds.set_data_matrix(dm_length)

		print msds.get_data_matrix()

	if args.peptide_normalize:
		#kdrew: input error checking
		assert (args.peptide_digest_format in ["raw", "counts"] )

		peptide_file = open(args.peptide_digest_filename, "rb") 
		peptide_counts = None
		if args.peptide_digest_format == "raw":
			peptide_counts = False
		if args.peptide_digest_format == "counts":
			peptide_counts = True
		dm_peptide = nu.normalize_peptide_count(msds.get_data_matrix(), msds.get_id_dict(), peptide_file, peptide_counts=peptide_counts)
		msds.set_data_matrix(dm_peptide)


	if args.threshold_normalize:
		msds_dm = msds.get_data_matrix()
		msds_dm = nu.threshold(msds_dm, threshold=args.threshold)
		msds.set_data_matrix(msds_dm)


	#print msds.get_id_dict()

	cPickle.dump( msds, open(args.out_filename,"wb"))

if __name__ == "__main__":
	main()

