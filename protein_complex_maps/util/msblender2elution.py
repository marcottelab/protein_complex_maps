

import argparse
import numpy as np
import pandas as pd
import os.path

def main():

	parser = argparse.ArgumentParser(description="Tool to read in msblender file and produce elution file")
	parser.add_argument("--prot_count_files", action="store", dest="prot_count_files", nargs='+', required=True, 
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--add_total_count", action="store_true", dest="add_total_count", required=False, default=False,
						help="Original file format has the first column as the sum of all the columns, default=False")
	parser.add_argument("--spectral_count_type", action="store", dest="spectral_count_type", required=False, default='Unique',
						help="Which spectral count type should be stored? Unweighted, Weighted or Unique (default)")

	args = parser.parse_args()

	elution_profile_df = pd.DataFrame() 
	for filename in args.prot_count_files:
		df = pd.read_table(filename, index_col=0)

		fraction_name = os.path.basename(filename).split('.')[0]
		print fraction_name
		single_elution_df = pd.DataFrame(df[args.spectral_count_type])
		single_elution_df.columns = [fraction_name]
		print single_elution_df
		elution_profile_df = pd.concat([elution_profile_df, single_elution_df], axis=1)


	elution_profile_df = elution_profile_df.fillna(0)

	if args.add_total_count:
		total_count = pd.DataFrame(elution_profile_df.sum(axis=1))
		total_count.columns = ['TotalCount']
		elution_profile_df = pd.concat([total_count, elution_profile_df], axis=1)

	elution_profile_df.to_csv(args.output_filename, sep='\t')



if __name__ == "__main__":
	main()



