
import sys
import numpy as np
import argparse
import pickle
import pandas
import gc
#from guppy import hpy

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu


def main():

	#hmem = hpy()

	parser = argparse.ArgumentParser(description="Computes correlation of input msds and outputs to file")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, default=None, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
						help="Ignore missing protein ids in msds")
	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, default=None,
						help="Threshold of correlation coefficients to store (ex. 0.0 = all corrcoef greater than 0.0), default=None")
	parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
						help="Filename of output file, default=None which prints to stdout")
	parser.add_argument("--remove_self", action="store_true", dest="remove_self", required=False, default=False, 
						help="Remove results where correlation is between the same indices")
	parser.add_argument("--offset", action="store", type=int, dest="offset", required=False, default=0, 
						help="""Add offset to indices, (ex. apcluster requires 1-indexing, so use --offset=1), 
								WARNING: don't forget to convert back when done, default=0""")
	parser.add_argument("--pieces", action="store", type=int, dest="pieces", required=False, default=10, 
						help="""To save memory break data into pieces, default=10""")

	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )
	print "after msds"
	#print hmem.heap()

	#if args.proteins != None:
	#	data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 
	#else:
	#	data_set = msds.get_data_matrix()
	#	new_id_map = msds.get_name2index()

	#kdrew: get data frames
	dataframe = msds.get_data_frame()
	print "after dataframe"
	#print hmem.heap()

	del msds
	gc.collect()
	print "del msds"
	#print hmem.heap()

	#dataframe = pandas.concat([dataframe,df])

	#kdrew: run correlation
	corrMat = dataframe.corr()
	print "after corrMat"
	#print hmem.heap()

	del dataframe
	gc.collect()
	print "del dataframe"
	#print hmem.heap()

	#kdrew: make indices
	column_list = corrMat.columns.tolist() 
	index_list = corrMat.index.tolist()
	#kdrew: make sure columns and indices are the same
	assert column_list == index_list
	index_number_list = range( args.offset, len(index_list)+args.offset )

	index_map = [(index_number_list[i],name) for i, name in enumerate(index_list)]

	corrMat.columns = index_number_list
	corrMat.index = index_number_list

	corrMat_stack = corrMat.stack()
	print "after corrMat_stack"
	#print hmem.heap()

	del corrMat
	gc.collect()
	print "del corrMat"
	#print hmem.heap()

	if args.remove_self:
		#kdrew: use generator instead of list
		corrMat_stack = corrMat_stack.drop(( (k,k) for k in index_list ) )
	print "after removed self"

	##kdrew: break stack into pieces so hopefully save some memory
	#corrMat_pieces = np.array_split(corrMat_stack, args.pieces)
	#new_corrMat_pieces = []

	#for stack_piece in corrMat_pieces:
	#	#kdrew: remove entries where both the column and index are the same
	#	if args.remove_self:
	#		#corrMat_stack = corrMat_stack[[ k for k in corrMat_stack.index if k[0] != k[1]]]
	#		#stack_piece = stack_piece.drop([ (k,k) for k in index_list ] )
	#		#stack_piece = stack_piece[[ k for k in stack_piece.index if k[0] != k[1]]]
	#		stack_piece = stack_piece[( k for k in stack_piece.index if k[0] != k[1] )]
	#	
	#	print "after removed self"
	#	#print hmem.heap()

	#	#kdrew: threshold entries by input
	#	#corrMat_stack = corrMat_stack[corrMat_stack >= args.threshold]
	#	if args.threshold != None:
	#		stack_piece = stack_piece[stack_piece >= args.threshold]
	#		print "after threshold entries"
	#	
	#	new_corrMat_pieces.append(stack_piece)
	#	gc.collect()
    #
	#
	#corrMat_stack = pandas.concat(new_corrMat_pieces)
	
	#print "after threshold entries"
	#print hmem.heap()

	out_handle = None
	if args.out_filename == None:
		out_handle = sys.stdout
		map_out_handle = sys.stdout
	else:
		out_handle = open(args.out_filename,"wb")
		map_out_handle = open(args.out_filename+".map","wb")

	#for i,j in results:
	#	if i == j and args.remove_self:
	#		continue
	#	out_handle.write("%s\t%s\t%s\n" % (i+args.offset, j+args.offset, results[(i,j)],))

	corrMat_stack.to_csv(out_handle, sep=' ')
	print "after to_csv"
	#print hmem.heap()

	for i,j in index_map:
		map_out_handle.write("%s %s\n" % (i,j))

	out_handle.close()

if __name__ == "__main__":
	main()


