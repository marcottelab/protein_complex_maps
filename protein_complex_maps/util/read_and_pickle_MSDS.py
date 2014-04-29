

import glob
import cPickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.peptide_util as ppu


def main():

	parser = argparse.ArgumentParser(description="Tool to read in Mass Spec Data Set (MSDS) from elution files and pickle it")
	parser.add_argument("--input_msds_files", action="store", dest="msds_filenames", nargs='+', required=False, default=None,
						help="Filenames of elution files")
	parser.add_argument("--input_peplist_files", action="store", dest="peplist_filenames", nargs='+', required=False, default=None,
						help="Filenames of peptide list files")
	parser.add_argument("--input_pepdict_file", action="store", dest="pepdict_filename", required=False, default=None,
						help="Filename of peptide list files")
	parser.add_argument("--protein_list", action="store", dest="protein_list", nargs='+', required=False, default=None,
						help="List of protein identifiers")
	parser.add_argument("--protein_list_file", action="store", dest="protein_list_file", required=False, default=None,
						help="File of protein identifiers")
	parser.add_argument("--spectral_count_mean", action="store_true", dest="spectral_count_mean", required=False, default=False,
						help="Set the spectral count to the mean of all peptide spectral counts for a protein")
	parser.add_argument("--peptide_threshold", action="store", type=int, dest="peptide_threshold", required=False, default=1,
						help="Protein must have specified number of peptides identified or else value is set to 0.0")
	parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
						help="Output filename of MSDS pickle")

	args = parser.parse_args()

	#kdrew: read in ms files
	#sample_filenames = "/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_*"

	#kdrew: error checking
	if not args.msds_filenames:
		if not args.peplist_filenames or not args.pepdict_filename:
			print "\nError: Specify either --input_msds_files or (--input_peplist_files and --input_pepdict_file)\n"
			parser.print_help()
			return
		

	msds = rd.MSDataSet()

	#kdrew: add to msds by elution files
	if args.msds_filenames:
		#for sample_filename1 in glob.iglob(sample_filenames):
		for sample_filename1 in args.msds_filenames:
			print sample_filename1
			sample_file1 = open(sample_filename1, 'rb')
			msds.load_file(sample_file1, header=True)

			#data_matrix = msds.get_data_matrix()

	#kdrew: add to msds by peptide counts
	if args.peplist_filenames and args.pepdict_filename:
		pepdict_file = open(args.pepdict_filename,'rb')

		if args.protein_list:
			peptide_dict = ppu.read_peptide_dict(pepdict_file, args.protein_list)
		elif args.protein_list_file:
			lfile = open(args.protein_list_file)
			protein_list = []
			for line in lfile.readlines():
				protein_list.append(line.split()[0])
			peptide_dict = ppu.read_peptide_dict(pepdict_file, protein_list)
		else:
			peptide_dict = ppu.read_peptide_dict(pepdict_file)

		for peptide_list_filename in args.peplist_filenames:
			print peptide_list_filename
			peptide_list_file = open(peptide_list_filename, 'rb')
			peptide_count_dict = ppu.read_peptide_list(peptide_list_file)

			#kdrew: collate peptide counts with peptide protein mapping
			protein_count_dict = ppu.protein_counts_by_peptide(peptide_dict, peptide_count_dict)

			#kdrew: add protein spectral counts to msds
			msds.create_by_peptide_counts( protein_count_dict, args.spectral_count_mean, int(args.peptide_threshold))

	print msds.get_id_dict()

	cPickle.dump( msds, open( args.out_filename, "wb"))

if __name__ == "__main__":
	main()




