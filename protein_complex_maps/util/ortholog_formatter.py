

import glob
import cPickle
import argparse


def main():

	parser = argparse.ArgumentParser(description="Tool to read in Mass Spec Data Set (MSDS) from elution files and pickle it")
	parser.add_argument("--filename", action="store", dest="filename", required=True,
						help="Filenames of Blake's inparanoid files")
	parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
						help="Output filename ")
	parser.add_argument("--species1", action="store", dest="species1", required=True, 
						help="species listed in OrtoA column")
	parser.add_argument("--species2", action="store", dest="species2", required=True, 
						help="species listed in OrtoB column")

	args = parser.parse_args()

	#kdrew: read in ms files
	#sample_filenames = "/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_*"

	#kdrew: error checking
	if not args.filename:
		print "\nError: Specify --filename \n"
		parser.print_help()
		return

	outfile = open(args.out_filename, "wb")
		
	file1 = open(args.filename, 'rb')
	for line in file1.readlines():
		col_split = line.split('\t')
		cluster_id = col_split[0]
		guid = col_split[1]
		species1_list = col_split[2].split()
		species2_list = col_split[3].split()

		for prot, score in zip(species1_list[0::2], species1_list[1::2]):
			outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cluster_id, guid, args.species1, score, prot,))

		for prot, score in zip(species2_list[0::2], species2_list[1::2]):
			outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cluster_id, guid, args.species2, score, prot,))


	file1.close()
	outfile.close()

if __name__ == "__main__":
	main()




