

import numpy as np

import argparse
import pickle
import itertools as it
import scipy.io
import pandas as pd

import protein_complex_maps.protein_util as pu

def main():

	parser = argparse.ArgumentParser(description="Map network from id type to another id type and output sparse format (id1 id2 score)")
	parser.add_argument("--matrix", action="store", dest="matrix", required=True,
						help="Filename of matrix")
	parser.add_argument("--matrix_ids", action="store", dest="matrix_ids", required=True, 
						help="Filename of ids in matrix")
	parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
						help="Filename of output file")
	parser.add_argument("--map_pickle", action="store_true", dest="map_pickle", required=False, default=False,
						help="Flag to load map from pickle file")
	parser.add_argument("--map_pickle_file", action="store", dest="map_pickle_file", required=False, default="map.p",
						help="Filename of pickle map")
	parser.add_argument("--map_id_from", action="store", dest="map_id_from", required=False, default="ENSEMBL_ID",
						help="Map ids of this type to another type, default=ENSEMBL_ID (list can be seen http://www.uniprot.org/faq/28)")
	parser.add_argument("--map_id_to", action="store", dest="map_id_to", required=False, default="ACC",
						help="Map ids to this type, default=ACC (list can be seen http://www.uniprot.org/faq/28)")
	parser.add_argument("--no_map", action="store_true", dest="no_map", required=False, default=False,
						help="Use matrix ids as identifer")
	parser.add_argument("--input_matrix_format", action="store", dest="input_matrix_format", required=False, default="dense",
						help="Format of input matrix file, dense, dataframe or mm, default=dense")

	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, default=0.0,
						help="Only include edges that pass a min threshold (>=), default = 0.0")
	args = parser.parse_args()


	mb_set = set()

        matrix_id_list = None
    
        if args.input_matrix_format == "mm":
            mat = scipy.io.mmread(args.matrix).todense()
        elif args.input_matrix_format == "dense":
            mat = np.loadtxt(args.matrix)
        elif args.input_matrix_format == "dataframe":
            mat = pd.read_csv(args.matrix,index_col=0, sep='\t')
            matrix_id_list = mat.index
        else:
            print "incorrect matrix format: %s" % args.input_matrix_format
            return -1

        if matrix_id_list == None:
            matrix_id_list = []

            #kdrew: read in ids
            matrix_ids_filename = args.matrix_ids
            matrix_id_file = open(matrix_ids_filename,"rb")
            for line in matrix_id_file.readlines():
                line_id = line.strip()
                matrix_id_list.append(line_id)
            matrix_id_file.close()

        if not args.map_pickle and not args.no_map:
            #kdrew: use uniprot to map ids
            matrix_id2ACC_map = pu.map_protein_ids(matrix_id_list, args.map_id_from, "ACC")
            flatten_list = [item for sublist in matrix_id2ACC_map.values() for item in sublist]
            matrix_ACC2id_map = pu.map_protein_ids(flatten_list, "ACC", args.map_id_to)
            #kdrew: reconcile maps

            matrix_id_map = dict()
            #kdrew: for every inputted id
            for id1 in matrix_id2ACC_map: 
                #kdrew: get all uniprot ids
                for acc in matrix_id2ACC_map[id1]:
                    #kdrew: for all uniprot ids get translated id
                    for id2 in matrix_ACC2id_map[acc]:
                        #kdrew: if there is an in_id -> translated_id, and translated_id is not already used, store result in map
                        if id2 != None and id2 not in matrix_id_map.values():
                            if matrix_id_map.has_key(id1):
                                #kdrew: ignore if there are duplicate translated_ids, defaults to the first one
                                print "duplicate for %s: %s" % (id1, id2)
                            else:
                                matrix_id_map[id1] = id2
                        #kdrew: if no result or translated_id already taken, store original id in map
                        else:
                            matrix_id_map[id1] = id1

            pickle.dump(matrix_id_map, open(args.map_pickle_file,"wb"))
        elif args.no_map:
            matrix_id_map = dict()
            for id1 in matrix_id_list:
                matrix_id_map[id1] = id1
        else:
            matrix_id_map = pickle.load(open(args.map_pickle_file,"rb"))


        print mat
        index=[matrix_id_map[x] if matrix_id_map.has_key(x) else x for x in matrix_id_list ]
        #columns=[matrix_id_map[x] if matrix_id_map.has_key(x) else x for x in matrix_id_list ]
        df = pd.DataFrame(mat)
        df.index = index
        df.columns = index
        print df
        #print len(df.index)
        #print len(df.columns)
        #print '11218' in df.columns

        outfile = open(args.output_filename, "wb")
        for i,j in it.combinations(df.index, 2):
            #print i,j
            #print df.ix[i,j]
            if df.ix[i,j] >= args.threshold:
                outfile.write("%s\t%s\t%s\n" % (i,j,df.ix[i,j]))


if __name__ == "__main__":
	main()




