

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
	parser.add_argument("--add_total_count", action="store_true", dest="add_total_count", required=False, default=False,
						help="Original file format has the first column as the sum of all the columns, default=False")
	parser.add_argument("--spectral_count_type", action="store", dest="spectral_count_type", required=False, default='Unique',
						help="Which spectral count type should be stored? Unweighted, Weighted or Unique (default)")
	parser.add_argument("--group", action="store_true", dest="group", required=False, default=False,
						help="Flag for reading in msblender .group files, default=False")

	args = parser.parse_args()

	elution_profile_df = pd.DataFrame() 
	all_group_id_list = []
	for filename in args.prot_count_files:
		fraction_name = os.path.basename(filename).split('.')[0]
		print fraction_name

		if args.group:
			df = pd.read_table(filename, index_col=0)
		else:
			#kdrew: the non-group files have an extra line at the bottom for total counts
			df = pd.read_table(filename, index_col=0, skip_footer=1)

		if args.group:
			group = df['ProtID in Group']
			#kdrew: hacky bullshit to condense the index into the same list as the other protein ids in group
			values = group.fillna('').values.tolist()
			values = [x.split(';;') for x in values]
			print group.index
			group_id_list = zip(*[group.index,values])
			#kdrew: filter removes the empty strings added above from the fillna, the [[]]+[] adds string to list
			group_id_list = [ filter( None, [l[0]]+l[1] ) for l in group_id_list ]  
			print "flatten"
			print group_id_list
			all_group_id_list = all_group_id_list + group_id_list

		single_elution_df = pd.DataFrame(df[args.spectral_count_type])
		single_elution_df.columns = [fraction_name]
		#print single_elution_df
		elution_profile_df = pd.concat([elution_profile_df, single_elution_df], axis=1)


	elution_profile_df = elution_profile_df.fillna(0)

	if args.add_total_count:
		total_count = pd.DataFrame(elution_profile_df.sum(axis=1))
		total_count.columns = ['TotalCount']
		elution_profile_df = pd.concat([total_count, elution_profile_df], axis=1)

	elution_profile_df.to_csv(args.output_filename, sep='\t')

	if args.group:
		print "all_group_id_list"
		print len(all_group_id_list)
		group_id_set = set([tuple(x) for x in all_group_id_list])
		print len(group_id_set)
		group_id_set_dict = dict()
		for group_ids in group_id_set:
			for group_id in group_ids:
				for group_ids2 in group_id_set:
					if group_id in group_ids2:
						try:
							group_id_set_dict[group_id].add(group_ids2)
						except KeyError:
							#kdrew: not sure why another set is necessary, I thought the set above would have dealt with any redundancies
							group_id_set_dict[group_id] = set([group_ids2])
		
						
						
		print "group_id_set_dict"
		print group_id_set_dict
		print "length of sets"
		print Counter([len(group_id_set_dict[x]) for x in group_id_set_dict])


if __name__ == "__main__":
	main()



