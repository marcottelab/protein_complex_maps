
import sys
sys.setrecursionlimit(10000)
import numpy as np
import argparse
import pickle
import pandas
import gc
#from guppy import hpy

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu
import protein_complex_maps.normalization_util as nu


def main():

	parser = argparse.ArgumentParser(description="Outputs fractionation matrix to file")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, default=None, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
						help="Ignore missing protein ids in msds")
	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, default=0.0,
						help="Threshold of spectral counts to store, default=0.0")
	parser.add_argument("--binary", action="store_true", dest="binary", required=False, default=False,
                                                help="Flag to output binary: 1 if present (>= threshold) and 0 if absent (< threshold)")
	parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
						help="Filename of output file, default=None which prints to stdout")
	parser.add_argument("--dataframe", action="store_true", dest="dataframe", required=False, default=False, 
						help="Flag to output dataframe with row and column names")
	parser.add_argument("--break2blocks", action="store_true", dest="break2blocks", required=False, default=False, 
						help="Break up dataset into blocks based on correlation above defined threshold, default = False")
	parser.add_argument("--blocks_threshold", action="store", type=float, dest="blocks_threshold", required=False, default=0.5, 
						help="Correlation threshold to break dataset into blocks, default = 0.5")
	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )
	print "after msds"



	if args.proteins != None:
            data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 
	else:
            data_set = msds.get_data_matrix()
            new_id_map = msds.get_name2index()

        if args.break2blocks:
            #kdrew: create correlation matrix
            #cormat = squareform(pdist(data_set, "correlation"))
            cormat = np.corrcoef(data_set)
            #kdrew: cluster
            clusters = []
            #kdrew: create list of ids to manipulate
            ids = new_id_map.keys()[:]
            while len(ids) > 0:
                #kdrew: generate new data_sets
                cur_i = ids.pop()
                #print "cur_i: %s" % cur_i
                #print "cormat array: %s" % cormat[cur_i]
                cluster = correlation_indices(cur_i, cormat, args.blocks_threshold, set([cur_i]))
                print "cluster: %s" % cluster 
                clusters.append(cluster)
                for j in cluster:
                    try:
                        ids.remove(j)
                    except ValueError:
                        continue

            for cluster in clusters:
                proteins = []
                for i in cluster:
                    proteins.append(new_id_map[i])
                print proteins

                



        if args.binary:
            if args.dataframe:
                bin_data_set = nu.binary(data_set, args.threshold)
                msds.set_data_matrix(bin_data_set)
                msds.get_data_frame().to_csv(args.out_filename)
            else:
                data_set = nu.binary(data_set, args.threshold)
                np.savetxt(args.out_filename, data_set, fmt='%d')
        else:
            if args.dataframe:
                msds.get_data_frame().to_csv(args.out_filename)
            else:
                np.savetxt(args.out_filename, data_set)

        map_out_handle = open(args.out_filename+".map","wb")
        
        for i in range(len(new_id_map)):
            map_out_handle.write("%s\n" % (new_id_map[i]))
        
	map_out_handle.close()

def correlation_indices(i, cormat, threshold, processed=set()):
    #print "correlation_index i: %s" % i
    cor_ids = np.nonzero( cormat[i] > threshold )
    #print "cor_ids: %s" % cor_ids
    cor_ids_return_set= set(cor_ids[0])
    processed.add(i)
    if len(cor_ids[0]) > 0:
        for j in cor_ids[0]:
            if j not in processed:
                processed.add(j)
                #print "j:%s" % j
                cor_ids_return_set = cor_ids_return_set.union(correlation_indices(j,cormat,threshold, processed))

    cor_ids_return_set.add(i)
    return cor_ids_return_set 
    
    

if __name__ == "__main__":
	main()


