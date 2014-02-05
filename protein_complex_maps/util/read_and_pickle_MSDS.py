

import glob
import cPickle
import argparse

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu


def main():

	parser = argparse.ArgumentParser(description="Tool to read in Mass Spec Data Set (MSDS) from elution files and pickle it")
	parser.add_argument("--input_msds_files", action="store", dest="msds_filenames", nargs='+', required=True, 
						help="Filenames of elution files")
	parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
						help="Output filename of MSDS pickle")

	args = parser.parse_args()

	#kdrew: read in ms files
	#sample_filenames = "/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_*"

	msds = rd.MSDataSet()
	#for sample_filename1 in glob.iglob(sample_filenames):
	for sample_filename1 in args.msds_filenames:
		print sample_filename1
		sample_file1 = open(sample_filename1, 'rb')
		msds.load_file(sample_file1, header=True)

		data_matrix = msds.get_data_matrix()

	print msds.get_id_dict()

	cPickle.dump( msds, open( args.out_filename, "wb"))

if __name__ == "__main__":
	main()

