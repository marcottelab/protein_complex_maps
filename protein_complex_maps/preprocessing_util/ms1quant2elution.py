

import argparse
import numpy as np
import pandas as pd
import os.path
from collections import Counter

def main():

	parser = argparse.ArgumentParser(description="Tool to read in msblender file and produce elution file")
	parser.add_argument("--prot_count_files", action="store", dest="prot_count_files", nargs='+', required=True, 
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--fraction_name_from_filename", action="store_true", dest="fraction_name_from_filename", required=False, default=False,
						help="Flag for parsing fraction name from filename, default=False, defaults to sequential numbering")
	parser.add_argument("--parse_uniprot_id", action="store_true", dest="parse_uniprot_id", required=False, default=False,
						help="Parse out the uniprot id from the index and store in index, default=False")
	parser.add_argument("--msblender_format", action="store_true", dest="msblender_format", required=False, default=False,
						help="Original file format has the first column as the sum of all the columns, default=False")

	args = parser.parse_args()

	elution_profile_df = pd.DataFrame() 
	for i, filename in enumerate(args.prot_count_files):
            df = pd.read_csv(filename, index_col=0)

            single_elution_df = pd.DataFrame(df['abundance'])


            if args.fraction_name_from_filename:
                fraction_name = os.path.basename(filename).split('.')[0]
                print fraction_name

                single_elution_df.columns = [fraction_name]
            else:
                single_elution_df.columns = [i]

            #print single_elution_df
            elution_profile_df = pd.concat([elution_profile_df, single_elution_df], axis=1)


	elution_profile_df = elution_profile_df.fillna(0)


        if args.parse_uniprot_id:
            uniprots = [x.split('|')[1] if len(x.split('|'))>1 else x for x in elution_profile_df.index ]
            elution_profile_df.index = uniprots

	if args.msblender_format:
            total_count = pd.DataFrame(elution_profile_df.sum(axis=1))
            total_count.columns = ['TotalCount']
            elution_profile_df = pd.concat([total_count, elution_profile_df], axis=1)
            elution_profile_df.index.name='#ProtID'
            print elution_profile_df.index.name

	print elution_profile_df.index.name
        elution_profile_df.to_csv(args.output_filename, sep='\t')


if __name__ == "__main__":
	main()



