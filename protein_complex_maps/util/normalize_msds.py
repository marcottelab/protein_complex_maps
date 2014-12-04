

import glob
import cPickle
import pickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.stoichiometry.stoichiometry as st
import protein_complex_maps.stoichiometry.relative_stoichiometry as rs
import protein_complex_maps.protein_util as pu


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

	parser.add_argument("--peptide_digest_delimiter", action="store", dest="delimiter", required=False, default='|',
						help="Delimiter of peptide digest file, seperates multiple proteins")

	parser.add_argument("--threshold_normalize", action="store_true", dest="threshold_normalize", required=False, default=False,
						help="Normalize by threshold")

	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, 
						help="Threshold for which to normalize by")

	parser.add_argument("--standardize", action="store_true", dest="standardize", required=False, default=False,
						help="Standardize each data matrix (x - x_mean) / x_std ")

	parser.add_argument("--map_ids", action="store_true", dest="map_ids", required=False, default=False,
						help="Map one id type to another, set using map_id_from and map_id_to ")

	parser.add_argument("--map_id_from", action="store", dest="map_id_from", required=False, default="ENSEMBL_ID",
						help="Map ids of this type to another type, default=ENSEMBL_ID (list can be seen http://www.uniprot.org/faq/28)")

	parser.add_argument("--map_id_to", action="store", dest="map_id_to", required=False, default="ACC",
						help="Map ids to this type, default=ACC (list can be seen http://www.uniprot.org/faq/28)")

	parser.add_argument("--transfer_map", action="store_true", dest="transfer_map", required=False, default=False,
						help="Transfer protein id map from one msds to another")

	parser.add_argument("--transfer_msds_pickle", action="store", dest="transfer_msds_filename", required=False, 
						help="Filename of MSDS pickle to be transfered")

	parser.add_argument("--map_by_genename", action="store_true", dest="map_by_genename", required=False, default=False,
						help="Map genename to ACC, can set organism")

	parser.add_argument("--organism", action="store", dest="organism", required=False, default="",
						help="Map genenames from specified organism")

	parser.add_argument("--map_orthologs", action="store_true", dest="map_orthologs", required=False, default=False,
						help="Map orthologs from species1 to species2")
	parser.add_argument("--species1", action="store", dest="species1", required=False, default=None, 
						help="species listed in msds")
	parser.add_argument("--species2", action="store", dest="species2", required=False, default=None,
						help="species to map orthologs to")

	args = parser.parse_args()

	#length_normalize = False
	#threshold_normalize = False

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )
	if args.map_ids:
		msds.map_ids(args.map_id_from, args.map_id_to)

	if args.map_by_genename:
		msds.map_ids_by_genename(args.organism)

	if args.transfer_map:
		msds2 = pickle.load( open( args.transfer_msds_filename, "rb" ) )
		msds.transfer_map(msds2)

	if args.map_orthologs:
		ortholog_map = pu.get_ortholog( msds.get_name_list(), args.species1, args.species2, version="_origid", database="inparanoid_blake")
		msds.add_mapping( ortholog_map, args.species2+"_ortho" )

	#kdrew: normalize by length as well as map ids
	if args.length_normalize:
		print msds.get_id_dict()
		#kdrew: map ids
		msds.map_ids(args.map_id_from, args.map_id_to)
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
		dm_peptide = nu.normalize_peptide_count(msds.get_data_matrix(), msds.get_id_dict(), peptide_file, peptide_counts=peptide_counts, delimiter=args.delimiter)
		msds.set_data_matrix(dm_peptide)


	if args.threshold_normalize:
		msds_dm = msds.get_data_matrix()
		msds_dm = nu.threshold(msds_dm, threshold=args.threshold)
		msds.set_data_matrix(msds_dm)


	if args.standardize:
		msds_dm = msds.get_data_matrix()
		msds_dm = nu.standardize(msds_dm)
		msds.set_data_matrix(msds_dm)

	#print msds.get_id_dict()

	pickle.dump( msds, open(args.out_filename,"wb"))

if __name__ == "__main__":
	main()

