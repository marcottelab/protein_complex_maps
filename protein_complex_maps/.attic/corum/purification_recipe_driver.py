
import csv
import pickle
import argparse
import os.path
import sys
import logging

import protein_complex_maps.purification_recipe.purification_recipe as pr

def main():

	parser = argparse.ArgumentParser(description="Driver script for creating purification scripts for corum complexes")
	parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True,
							help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--complex_file", action="store", dest="complex_filename", required=False, default=None,
							help="Filename of complex protein identifiers ie. corum csv")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
							help="Filename of output plot")
	parser.add_argument("--fractionation_type", action="store", dest="fractionation_type_file", required=False, default="./fractionation_type.txt",
							help="File that describes the type of column used for a fractionation experiment")
	parser.add_argument("--protein_percent_threshold", action="store", dest="protein_percent_threshold", type=float, required=False, default=0.9,
							help="Percentage of proteins that need to be present in input fractions expressed as decimal, default 0.9 ie 90%" )

	args = parser.parse_args()

	complexes_dict = {}
	complex_names_dict = {}

	if args.complex_filename:
		#kdrew: read in complexes from corum
		with open(args.complex_filename, 'rb') as csvfile:
			complexes = csv.reader(csvfile, delimiter=';')
			for row in complexes:
				if row[3] == "Human" and row[4] != '':
					ids = row[4].split(',')
					#kdrew: remove ( ) around extra ids
					complexes_dict[row[0]] = [x.translate(None, '()') for x in ids]
					complex_names_dict[row[0]] = row[1]

	#kdrew: for every complex create a purification recipe
	for complex1 in complexes_dict:

		try:
			plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(complex1)+os.path.splitext(args.plot_filename)[1]
			print plot_fname

			print complexes_dict[complex1]

			pr_obj = pr.PurificationRecipe( args.msds_filenames, complexes_dict[complex1], args.protein_percent_threshold, plot_fname, plot_title="Purification Recipe: %s" % complex_names_dict[complex1], fractionation_type_file=args.fractionation_type_file)
			pr_obj.create_recipes()
			pr_obj.show_results()
			print pr_obj.has_results()
			if args.plot_filename != None and pr_obj.has_results():
				pr_obj.plot_results()


        #kdrew: if there is a problem continue to next complex
		except:
			#e = sys.exc_info()[0]
			#print "complex: %s error: %s" % (complex1,e)
			logging.exception("complex: %s" % (complex1))
			continue


if __name__ == "__main__":
    main()
