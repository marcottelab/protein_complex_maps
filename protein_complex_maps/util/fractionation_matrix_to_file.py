
import sys
import numpy as np
import argparse
import pickle
import pandas
import gc
#from guppy import hpy

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
	parser.add_argument("--map_name", action="store", dest="map_name", required=False, default=None, 
						help="Output mapped ids instead of master list ids, default=None")
	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )
	print "after msds"

	if args.proteins != None:
		data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 
	else:
		data_set = msds.get_data_matrix()
		new_id_map = msds.get_name2index()



        if args.binary:
            if args.dataframe:
                bin_data_set = nu.binary(data_set, args.threshold)
                msds.set_data_matrix(bin_data_set)
                msds.get_data_frame(map_ids=args.use_map_ids).T.to_csv(args.out_filename)
            else:
                data_set = nu.binary(data_set, args.threshold)
                np.savetxt(args.out_filename, data_set, fmt='%d')
        else:
            if args.dataframe:
                msds.get_data_frame(map_name=args.map_name).T.to_csv(args.out_filename)
            else:
                np.savetxt(args.out_filename, data_set)

        map_out_handle = open(args.out_filename+".map","wb")
        
        for i in range(len(new_id_map)):
            map_out_handle.write("%s\n" % (new_id_map[i]))
        
	map_out_handle.close()

if __name__ == "__main__":
	main()


