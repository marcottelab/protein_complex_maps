

import numpy as np
from numpy.random import normal
import matplotlib as mpl
mpl.use('Agg')

import argparse
import pickle
import MySQLdb
import itertools as it
from scipy.stats.stats import pearsonr
from scipy.stats import gaussian_kde
import scipy.io
import Bio.PDB
import matplotlib.pyplot as plt

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

import protein_complex_maps.protein_util as pu
import protein_complex_maps.pdb_util as pdbu

def main():

	parser = argparse.ArgumentParser(description="Plot precision recall of protein interactions for MB Stars graph estimation and correlation")
	parser.add_argument("--pdb_list_filename", action="store", dest="pdb_list_filename", required=False, default=None,
						help="File of benchmark pdbids")
	parser.add_argument("--pdbid", action="store", dest="pdbid", required=False, default=None,
						help="PDB id")
	parser.add_argument("--base_species", action="store", dest="base_species", required=False, default="Hsapiens",
						help="Species to evaluate, default=Hsapiens")
	parser.add_argument("--pisa_dir", action="store", dest="pisa_dir", required=True, default=None,
						help="Directory of pisa files for calculating surface area")
	parser.add_argument("--matrix_ids", action="store", dest="matrix_ids", required=False, default=None,
						help="Filename of ids in matrix")
	parser.add_argument("--mb_matrix_ids", action="store", dest="mb_matrix_ids", required=False, default=None,
						help="Filename of ids in mb matrix")
	parser.add_argument("--mb_map_pickle", action="store_true", dest="mb_map_pickle", required=False, default=False,
						help="Flag to load mb map from pickle file")
	parser.add_argument("--mb_map_pickle_file", action="store", dest="mb_map_pickle_file", required=False, default="mb_map.p",
						help="Filename of pickle map")
	args = parser.parse_args()


	pdb_list = []
	if args.pdb_list_filename != None:
            #kdrew: read in pdb list
            pdb_list_file = open(args.pdb_list_filename,"rb")
            for line in pdb_list_file.readlines():
                pdb_list.append(line.strip())
	elif args.pdbid != None:
            pdb_list.append(args.pdbid)
        else:
            print "Use either --pdb_list_filename or --pdbid to specify PDB"
            return


	interaction_list = []
	mb_interaction_list = []
	correlation_list = []
	nd_list = []
        mb_list = []
        biogrid_list = []

        #kdrew: read in MB ids
        mb_matrix_ids_filename = args.mb_matrix_ids
        mb_matrix_id_file = open(mb_matrix_ids_filename,"rb")
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
        #print mb_matrix_uniprot_map

        mb_matrix_id_dict = dict()
        mb_matrix_id_lists = []
        for pid in mb_matrix_ensg_list:
            try:
                mb_matrix_id_lists.append(mb_matrix_uniprot_map[pid])
            except IndexError:
                mb_matrix_id_lists.append(None)

        matrix_ids_filename = args.matrix_ids

        #kdrew: read in matrix (ND) ids
        matrix_id_file = open(matrix_ids_filename,"rb")
        matrix_id_list = []
        for line in matrix_id_file.readlines():
                matrix_id_list.append(line.strip())
        matrix_id_file.close()

        protein_mapped = dict()
        protein_ortho = dict()
        protein_chain = dict()
        protein_accs = dict()
        combinations = dict()
	#kdrew: for every pdb
	for pdbid in pdb_list:
            print pdbid

            chain2acc = pu.get_pdb_protein_ids(pdbid)
            print "chain2acc: %s" % chain2acc
            print "# of pdb acc: %s" % (len(set(chain2acc.values())))
            protein_accs[pdbid] =  list(set(chain2acc.values()))

            acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=args.base_species, reversible=True)
            print "acc2base_species: %s" % acc2base_species
            print "# of pdb ortho: %s" % (1.0*len(set(acc2base_species.keys()))/2.0)
            #kdrew: remember these are both species and base species, divide by two
            protein_ortho[pdbid] = acc2base_species.keys()


            pdb_matrix_id_list = []
            for pid in matrix_id_list:
                if acc2base_species.has_key(pid):
                    pdb_matrix_id_list.append(pid)
                else:
                    continue

            pdb_mb_matrix_id_list = []
            for i, pid_list in enumerate(mb_matrix_id_lists):
                notfound_flag = True
                for pid in pid_list:
                    if acc2base_species.has_key(pid) and notfound_flag:
                        pdb_mb_matrix_id_list.append(pid)
                        mb_matrix_id_dict[pid] = i
                        notfound_flag = False
                        print "%s %s" % (pid, mb_matrix_id_dict[pid])

            protein_mapped[pdbid] = pdb_mb_matrix_id_list

            #kdrew: load pisa file
            pisaInt = pdbu.PISA_Interfaces( args.pisa_dir+'/'+pdbid+'_pisa', pdbid=pdbid )

            combinations[pdbid] = 0
            for acc1, acc2 in it.combinations(pdb_mb_matrix_id_list,2):
                surface_area = pisaInt.surface_area_by_acc(acc1, acc2, base_species=args.base_species)
                if surface_area != None:
                    interaction_list.append(1)
                else:
                    interaction_list.append(0)

                combinations[pdbid] += 1

                #nd_list.append( nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])
                #correlation_list.append( corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])
                #biogrid_list.append( biogrid_mat[biogrid_id_list.index(acc1), biogrid_id_list.index(acc2)])
                #mb_list.append( mb_mat[mb_matrix_id_dict[acc1], mb_matrix_id_dict[acc2]])


        for pdbid in protein_accs.keys():
            print pdbid
            print len(set(protein_accs[pdbid]))
            print len(set(protein_ortho[pdbid]))
            print len(set(protein_mapped[pdbid]))
            print "combinations: %s" % combinations[pdbid]

        print "total interactions: %s" % sum(interaction_list)

if __name__ == "__main__":
	main()




