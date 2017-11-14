import os.path
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

def main():

	parser = argparse.ArgumentParser(description="Tool to read in msblender file and produce elution file")
	parser.add_argument("--prot_count_files", action="store", dest="prot_count_files", nargs='+', required=True, 
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
						help="Filenames of MSblender prot_count files")
	parser.add_argument("--spectral_count_type", action="store", dest="spectral_count_type", required=False, default='Unique',
						help="Which spectral count type should be stored? Unweighted, Weighted or Unique (default)")
	parser.add_argument("--remove_zero_unique", action="store_true", dest="remove_zero_unique", required=False, default=False,
						help="Flag for removing entries with zero unique matches, default=False")
	parser.add_argument("--fraction_name_from_filename", action="store_true", dest="fraction_name_from_filename", required=False, default=False,
						help="Flag for parsing fraction name from filename, default=False, defaults to sequential numbering")
	
	# Flags to write a clustered and annotated outfile
	parser.add_argument("--write_verbose_file", action="store_true",dest="write_verbose_file",required=False, default=False,
						help="If True, also write an tsv file with row and column sums, clustered on rows. Can also add annotations with --path_to_annotations")
	parser.add_argument("--hclust_method", action="store",dest="hclust_method",required=False,default='average',choices=["single","complete","average","weighted","centroid","median","ward"],
						help="Method for hierarchical clustering of rows. Only activated if --write_verbose_file if True")
	parser.add_argument("--hclust_metric", action="store",dest="hclust_metric",required=False,default="euclidean",choices=["euclidean","canberra","pearson","spearman"],
						help="Distance or correlation metric used for clustering")
	parser.add_argument("--path_to_annotations", action='store', dest='path_to_annotations', required=False, default=None,
						help="Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second. Must also specify --write_verbose_file")
	
	# Formatting flags
	parser.add_argument("--contaminant_flag", action='store', dest='contaminant_flag', required=False, default="CONTAMINANT",
						help="If specified, remove rows with this string from output file. Can specify None if you want to retain all rows for some reason")
	parser.add_argument("--msblender_format", action="store_true", dest="msblender_format", required=False, default=False,
						help="Original file format has the first column as the sum of all the columns, default=False")
	parser.add_argument("--parse_uniprot_id", action="store_true", dest="parse_uniprot_id", required=False, default=False,
						help="Parse out the uniprot id from the index and store in index, default=False")

	args = parser.parse_args()
	
	##BJL: asserts
	if args.path_to_annotations != None and args.write_verbose_file == False:
		raise Exception("To add annotations, you must specify --write_verbose_file")
		assert os.path.isfile(args.path_to_annotations), "File {} not found".format(args.path_to_annotations)

	elution_profile_df = pd.DataFrame() 
	
	for i, filename in enumerate(args.prot_count_files):
		print "Reading {}".format(filename)
		df = pd.read_table(filename, index_col=0)
		
		## BJL: Useful code if we want to eventually parse other MSBlender outfiles besided .group
		#else:
		#	#kdrew: the non-group files have an extra line at the bottom for total counts
		#	df = pd.read_table(filename, index_col=0, skip_footer=1)
		
		single_elution_df = pd.DataFrame(df[args.spectral_count_type])

		if args.remove_zero_unique:
			single_elution_df = single_elution_df[single_elution_df['Unique'] > 0]

		if args.fraction_name_from_filename:
			fraction_name = os.path.basename(filename).split('.')[0]
			#print fraction_name # Debugging

			single_elution_df.columns = [fraction_name]
		else:
			single_elution_df.columns = [i]

		#print single_elution_df # Debugging
		elution_profile_df = pd.concat([elution_profile_df, single_elution_df], axis=1)


	elution_profile_df = elution_profile_df.fillna(0)
	
	if args.contaminant_flag != None:
		elution_profile_df = elution_profile_df[~elution_profile_df.index.str.contains(args.contaminant_flag)]

	if args.parse_uniprot_id:
		uniprots = [x.split('|')[1] if len(x.split('|'))>1 else x for x in elution_profile_df.index ]
		elution_profile_df.index = uniprots

	if args.msblender_format:
		total_count = pd.DataFrame(elution_profile_df.sum(axis=1))
		total_count.columns = ['TotalCount']
		elution_profile_df = pd.concat([total_count, elution_profile_df], axis=1)
		elution_profile_df.index.name='#ProtID'
		print elution_profile_df.index.name

	#print elution_profile_df.index.name # Debugging
	elution_profile_df.to_csv(args.output_filename, sep='\t')
	
	## Write an tab delimited file with row and column sums, with the rows clustered according to some metric
	## Can then add in sparklines for a nice viz of the elution data.
	if args.write_verbose_file: 
		
		print "Clustering method={}, metric={}".format(args.hclust_method,args.hclust_metric)
		
		# Unfortunately using correlations doesn't work well because the matrix correlation functions
		labels_index = pd.Series(elution_profile_df.index,index=range(len(elution_profile_df))) # Create indexer to manage scipy output
		
		if args.hclust_metric == "pearson":
			corr_mat = np.around( np.corrcoef(elution_profile_df), decimals = 4)
			dist_mat = (1. - corr_mat) / 2. # BJL: Invert correlations and move to 0-1 axis, i.e. convert them to distances
			distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
			link = hclust.linkage(distArray, method=args.hclust_method)
		elif args.hclust_metric == "spearman":
			corr_mat = np.around( stats.spearmanr(elution_profile_df.T)[0], decimals = 4 )
			dist_mat = (1. - corr_mat) / 2.
			distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
			link = hclust.linkage(distArray, method=args.hclust_method)
		elif args.hclust_metric in ("euclidean","canberra"):
			link = hclust.linkage(elution_profile_df, method=args.hclust_method, metric=args.hclust_metric)
		else:
			raise # shouldn't be necessary b/c arparse should catch it with "choices", but whatever
		
		leaves_list = hclust.leaves_list(link)
		ordered_index = labels_index.iloc[leaves_list] # reorder index from hier clustering
		ordered_index.to_csv("ordered_index.csv")
		
		if args.path_to_annotations != None:
			print "Adding annotations"
			annots = pd.read_table(args.path_to_annotations,index_col=0,header=None,names=["annotation"])
			#print "length of annotation df: {}".format(len(annots)) # Debugging
			annotated_df = annots.join(elution_profile_df,how='right')
		
		print "Summing"
		colSums = pd.Series(elution_profile_df.sum(),name="colSum") # add to df later
		annotated_df["rowSum"] = elution_profile_df.sum(axis=1)
		clustered_index_name = "_".join(["order",args.hclust_metric,args.hclust_method])
		annotated_df[clustered_index_name] = range(len(ordered_index))

		# reorder by hier clustering and ordering column
		clustered_df = annotated_df.reindex( ordered_index )
		clustered_df[clustered_index_name] = range(len(ordered_index))
		
		# reorder columns
		column_order = ["annotation","rowSum",clustered_index_name] + [i for i in clustered_df.columns if i not in ["annotation","rowSum",clustered_index_name]] # haaacky bleerggh (BJL)
		reordered_df = clustered_df[column_order]
		
		# add column sum row
		outdf = reordered_df.append(colSums) # BJL: Must add annots, then row sum, the cluster order, then col sums -- in that order
		
		print "Writing excel file"
		out_file = ".".join( args.output_filename.split(".")[:-1] + ["clustered","xlsx"] )
		outdf.to_excel(out_file)


if __name__ == "__main__":
	main()



