

import glob
import cPickle

import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu


def main():

	#kdrew: read in ms files
	sample_filenames = "/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_*"

	msds = rd.MSDataSet()
	for sample_filename1 in glob.iglob(sample_filenames):
		print sample_filename1
		sample_file1 = open(sample_filename1, 'rb')
		msds.load_file(sample_file1, header=True)

		data_matrix = msds.get_data_matrix()
		clean_data_matrix_normalized1 = nu.normalize_over_columns(data_matrix)
		msds.set_data_matrix(clean_data_matrix_normalized1)
	
	print msds.get_id_dict()

	cPickle.dump( msds, open( "HS_ms2_elutions_msds.p", "wb"))

if __name__ == "__main__":
	main()

