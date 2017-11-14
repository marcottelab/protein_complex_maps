

import numpy as np

import argparse
import pickle
import itertools as it
import scipy.io

import protein_complex_maps.protein_util as pu

def main():

	parser = argparse.ArgumentParser(description="Extract sub-network from full network given protein ids")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--genenames", action="store_true", dest="genenames", required=False, default=False,
						help="Set labels to genenames")
	parser.add_argument("--mb_matrix", action="store", dest="mb_matrix", required=False, default=None,
						help="Filename of mb matrix")
	parser.add_argument("--mb_matrix_ids", action="store", dest="mb_matrix_ids", required=False, default=None,
						help="Filename of ids in mb matrix")
	parser.add_argument("--corr_matrix", action="store", dest="corr_matrix", required=False, default=None,
						help="Filename of corr matrix")
	parser.add_argument("--matrix_ids", action="store", dest="matrix_ids", required=False, default=None,
						help="Filename of ids in matrix")
	parser.add_argument("--mb_map_pickle", action="store_true", dest="mb_map_pickle", required=False, default=False,
						help="Flag to load mb map from pickle file")
	parser.add_argument("--mb_map_pickle_file", action="store", dest="mb_map_pickle_file", required=False, default="mb_map.p",
						help="Filename of pickle map")
	parser.add_argument("--mb_uniprot", action="store_true", dest="mb_uniprot_flag", required=False, default=False,
						help="mb_matrix_ids has uniprot ids instead of ENSEMBL_IDs")
	parser.add_argument("--neighbors", action="store_true", dest="neighbors", required=False, default=False,
						help="Flag to grow out to neighboring nodes of input proteins")
	parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, default=0.0,
						help="Only include edges that pass a min threshold (>=), default = 0.0")
	parser.add_argument("--map_ids", action="store_true", dest="map_ids", required=False, default=False,
						help="Map one id type to another, set using map_id_from and map_id_to ")
	parser.add_argument("--map_id_from", action="store", dest="map_id_from", required=False, default="ENSEMBL_ID",
						help="Map ids of this type to another type, default=ENSEMBL_ID (list can be seen http://www.uniprot.org/faq/28)")
	parser.add_argument("--map_id_to", action="store", dest="map_id_to", required=False, default="ACC",
						help="Map ids to this type, default=ACC (list can be seen http://www.uniprot.org/faq/28)")
	args = parser.parse_args()


	mb_set = set()


	#kdrew: read in ids
	mb_matrix_ids_filename = args.mb_matrix_ids
	mb_matrix_id_file = open(mb_matrix_ids_filename,"rb")
	mb_matrix_id_lists = []
	mb_matrix_id_list = []

	if args.mb_uniprot_flag:
            for line in mb_matrix_id_file.readlines():
                #kdrew: needs to be in list because that is how it is handled below
                mb_matrix_id_lists.append([line.strip(),])
                mb_matrix_id_list.append(line.strip())

            mb_matrix_id_file.close()
	else:
            mb_matrix_ensg_list = []
            for line in mb_matrix_id_file.readlines():
                line_id = line.strip()
                mb_matrix_ensg_list.append(line_id)
            mb_matrix_id_file.close()
            if not args.mb_map_pickle:
                mb_matrix_uniprot_map = pu.map_protein_ids(mb_matrix_ensg_list, "ENSEMBL_ID", "ACC")
                pickle.dump(mb_matrix_uniprot_map, open(args.mb_map_pickle_file,"wb"))
            else:
                mb_matrix_uniprot_map = pickle.load(open(args.mb_map_pickle_file,"rb"))

            for pid in mb_matrix_ensg_list:
                try:
                    mb_matrix_id_lists.append(mb_matrix_uniprot_map[pid])
                except IndexError:
                    mb_matrix_id_lists.append(None)

                try:
                    #print mb_matrix_uniprot_map[pid][0]
                    mb_matrix_id_list.append(mb_matrix_uniprot_map[pid][0])
                except IndexError:
                    mb_matrix_id_list.append(None)

	mb_matrix_id_dict = dict()

        mb_mat = scipy.io.mmread(args.mb_matrix).todense()

        if args.corr_matrix != None:
            corr_mat = np.loadtxt(args.corr_matrix)

        if args.matrix_ids != None:
            matrix_id_file = open(args.matrix_ids,"rb")
            matrix_id_list = []
            for line in matrix_id_file.readlines():
                matrix_id_list.append(line.strip())
            matrix_id_file.close()

        for i, pid_list in enumerate(mb_matrix_id_lists):
            notfound_flag = True
            for pid in pid_list:
                if pid in args.proteins and notfound_flag:
                    mb_matrix_id_dict[pid] = i
                    notfound_flag = False

	if args.genenames:
            genename_map = pu.get_genenames_uniprot( args.proteins )
            print genename_map

	if args.map_ids:
		msds.map_ids(args.map_id_from, args.map_id_to)


        if args.neighbors:
            for acc1 in args.proteins:
                #kdrew: get all non-zero entries in mb_mat
                try:
                    whereobj = np.where(mb_mat[mb_matrix_id_dict[acc1],:] > 0.0)
                except KeyError:
                    continue
                for i in np.array(whereobj[1])[0]:
                    #kdrew: find acc2 
                    #print i
                    acc2 = mb_matrix_id_list[i]
                    #print "neighboring acc: %s to %s" % (acc2, acc1)
                    #kdrew: copied from above
                    for i, pid_list in enumerate(mb_matrix_id_lists):
                        notfound_flag = True
                        for pid in pid_list:
                            if pid == acc2 and notfound_flag:
                                mb_matrix_id_dict[pid] = i
                                #print pid
                                notfound_flag = False

                    

                    if acc2 != None:
                        if mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]] >= args.threshold:
                            mb_set.add((acc1,acc2,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]]))

            if args.genenames:
                #kdrew: get genenames for additional neighbor nodes
                genename_map.update(pu.get_genenames_uniprot( [i[1] for i in mb_set] ))

            for i in mb_set:
                if args.genenames:
                    try:
                        acc1_name = genename_map[i[0]]
                    except KeyError:
                        print "no genename for %s" % i[0]
                        acc1_name = i[0]
                    try:    
                        acc2_name = genename_map[i[1]]
                    except KeyError:
                        print "no genename for %s" % i[1]
                        acc2_name = i[1]
                print "%s %s %s" % (acc1_name, acc2_name,i[2])

        #kdrew: search only combinations of input proteins
        else:
            for acc1, acc2 in it.combinations(args.proteins,2):
                if args.genenames:
                    try:
                        acc1_name = genename_map[acc1]
                    except KeyError:
                        print "no genename for %s" % acc1
                        acc1_name = acc1
                    try:    
                        acc2_name = genename_map[acc2]
                    except KeyError:
                        print "no genename for %s" % acc2
                        acc2_name = acc2

                mb_set.add((acc1_name,acc2_name,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]],corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)]))



            for i in mb_set:
                print "%s %s %s %s" % i


if __name__ == "__main__":
	main()




