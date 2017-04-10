

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

import pandas as pd

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
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename of output plot")
    parser.add_argument("--pisa_dir", action="store", dest="pisa_dir", required=True, default=None,
                                            help="Directory of pisa files for calculating surface area")
    parser.add_argument("--nd_matrix", action="store", dest="nd_matrix", required=False, default=None,
                                            help="Filename of nd matrix")
    parser.add_argument("--mb_matrix", action="store", dest="mb_matrix", required=False, default=None,
                                            help="Filename of mb matrix")
    parser.add_argument("--corr_matrix", action="store", dest="corr_matrix", required=False, default=None,
                                            help="Filename of corr matrix")
    parser.add_argument("--matrix_ids", action="store", dest="matrix_ids", required=False, default=None,
                                            help="Filename of ids in matrix")
    parser.add_argument("--mb_matrix_ids", action="store", dest="mb_matrix_ids", required=False, default=None,
                                            help="Filename of ids in mb matrix")
    parser.add_argument("--biogrid_matrix", action="store", dest="biogrid_matrix", required=False, default=None,
                                            help="Filename of biogrid matrix")
    parser.add_argument("--biogrid_ids", action="store", dest="biogrid_ids", required=False, default=None,
                                            help="Filename of ids in matrix from biogrid")
    parser.add_argument("--replace_PDBID", action="store_true", dest="replace_PDBID", required=False, default=False,
                                            help="Flag to replace PDBID in nd_matrix, corr_matrix and matrix_ids args with pdb ids in pdb_list_filename")
    parser.add_argument("--mb_map_pickle", action="store_true", dest="mb_map_pickle", required=False, default=False,
                                            help="Flag to load mb map from pickle file")
    parser.add_argument("--mb_map_pickle_file", action="store", dest="mb_map_pickle_file", required=False, default="mb_map.p",
                                            help="Filename of pickle map")
    parser.add_argument("--mb_uniprot", action="store_true", dest="mb_uniprot_flag", required=False, default=False,
                                            help="mb_matrix_ids has uniprot ids instead of ENSEMBL_IDs")
    parser.add_argument("--corr_format", action="store", dest="corr_format", required=False, default="loadtxt",
                                            help="Format of the correlation matrix file, use loadtxt or read_csv to read in, default = loadtxt")
    args = parser.parse_args()

    print "start of script"

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

    pdb_mb_dict = {}

    #kdrew: read in ids
    mb_matrix_ids_filename = args.mb_matrix_ids
    mb_matrix_id_file = open(mb_matrix_ids_filename,"rb")
    mb_matrix_id_lists = []

    if args.mb_uniprot_flag:
        for line in mb_matrix_id_file.readlines():
            #kdrew: needs to be in list because that is how it is handled below
            mb_matrix_id_lists.append([line.strip(),])

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

        
        #kdrew: building reverse map from ACC to ENSEMBL_ID
        mb_matrix_uniprot_map_rev = dict()
        for key in mb_matrix_uniprot_map:
            for id1 in mb_matrix_uniprot_map[key]:
                if id1 != None:
                    mb_matrix_uniprot_map_rev[id1] = key

        #mb_matrix_uniprot_map = pu.map_protein_ids(mb_matrix_ensg_list, "ENSEMBL_ID", "ACC")
        #print mb_matrix_uniprot_map

        for pid in mb_matrix_ensg_list:
            try:
                mb_matrix_id_lists.append(mb_matrix_uniprot_map[pid])
            except IndexError:
                mb_matrix_id_lists.append(None)


    mb_matrix_id_dict = dict()

    corr_matrix_filename = args.corr_matrix
    if args.corr_format == "loadtxt":
        corr_mat = np.loadtxt(corr_matrix_filename)
    elif args.corr_format == "read_csv":
        #kdrew: indices are embedded in the file
        corr_mat = pd.read_csv(corr_matrix_filename, sep='\t', index_col=0, compression='gzip')
        print "Max: %s Min: %s" % (corr_mat.max(), corr_mat.min())
    else:
        print "unrecognized corr_format: %s" % args.corr_format
        return -1

    #kdrew: for every pdb
    for pdbid in pdb_list:
        print "for pdbid in pdb_list: %s"  % pdbid

        pdb_mb_dict[pdbid] = list()

        chain2acc = pu.get_pdb_protein_ids(pdbid)
        acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=args.base_species, reversible=True)

        matrix_ids_filename = args.matrix_ids
        biogrid_ids_filename = args.biogrid_ids
        nd_matrix_filename = args.nd_matrix
        mb_matrix_filename = args.mb_matrix
        biogrid_matrix_filename = args.biogrid_matrix

        if args.replace_PDBID:
            matrix_ids_filename = matrix_ids_filename.replace("PDBID",pdbid)
            biogrid_ids_filename = biogrid_ids_filename.replace("PDBID",pdbid)
            nd_matrix_filename = nd_matrix_filename.replace("PDBID",pdbid)
            mb_matrix_filename = mb_matrix_filename.replace("PDBID",pdbid)
            corr_matrix_filename = corr_matrix_filename.replace("PDBID",pdbid)
            biogrid_matrix_filename = biogrid_matrix_filename.replace("PDBID",pdbid)


        matrix_id_file = open(matrix_ids_filename,"rb")
        matrix_id_list = []
        for line in matrix_id_file.readlines():
            matrix_id_list.append(line.strip())
        matrix_id_file.close()

        #kdrew: read in biogrid ids
        if biogrid_ids_filename != None:
            biogrid_id_file = open(biogrid_ids_filename,"rb")
            biogrid_id_list = []
            for line in biogrid_id_file.readlines():
                biogrid_id_list.append(line.strip())
            biogrid_id_file.close()

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

        print pdb_mb_matrix_id_list
        print mb_matrix_id_dict

        #kdrew: readin matrices
        if nd_matrix_filename != None:
            nd_mat = np.loadtxt(nd_matrix_filename)
        if mb_matrix_filename != None:
            mb_mat = scipy.io.mmread(mb_matrix_filename).todense()


        if biogrid_matrix_filename != None:
            biogrid_mat = np.loadtxt(biogrid_matrix_filename)


        #kdrew: load pisa file
        pisaInt = pdbu.PISA_Interfaces( args.pisa_dir+'/'+pdbid+'_pisa', pdbid=pdbid )

        print "##############"
        #for acc1, acc2 in it.combinations(pdb_mb_matrix_id_list,2):
        for acc1, acc2 in it.combinations(set(pdb_mb_matrix_id_list),2):
            surface_area = pisaInt.surface_area_by_acc(acc1, acc2, base_species=args.base_species)
            print "acc1: %s acc2: %s sarea: %s" % (acc1, acc2, surface_area)
            if surface_area != None:
                interaction_list.append(1)
            else:
                interaction_list.append(0)

            if nd_matrix_filename != None:
                nd_list.append( nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])

            if args.corr_format == "loadtxt":
                correlation_list.append( corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])
            elif args.corr_format == "read_csv":
                correlation_list.append( corr_mat[mb_matrix_uniprot_map_rev[acc1]][mb_matrix_uniprot_map_rev[acc2]] )

            if biogrid_ids_filename != None:
                biogrid_list.append( biogrid_mat[biogrid_id_list.index(acc1), biogrid_id_list.index(acc2)])

        #for acc1, acc2 in it.combinations(pdb_mb_matrix_id_list,2):
        #    surface_area = pisaInt.surface_area_by_acc(acc1, acc2, base_species=args.base_species)

            is_interaction = None
            if surface_area != None:
                mb_interaction_list.append(1)
                is_interaction = 1
            else:
                mb_interaction_list.append(0)
                is_interaction = 0

            if mb_matrix_filename != None and nd_matrix_filename != None:
                mb_list.append( mb_mat[mb_matrix_id_dict[acc1], mb_matrix_id_dict[acc2]])
                if args.corr_format == "loadtxt":
                    pdb_mb_dict[pdbid].append((acc1,acc2,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]],is_interaction,surface_area,pdbid,nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)],corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)]))
                elif args.corr_format == "read_csv":
                    pdb_mb_dict[pdbid].append((acc1,acc2,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]],is_interaction,surface_area,pdbid,nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)],corr_mat[mb_matrix_uniprot_map_rev[acc1]][mb_matrix_uniprot_map_rev[acc2]]))

            elif mb_matrix_filename != None:
                mb_list.append( mb_mat[mb_matrix_id_dict[acc1], mb_matrix_id_dict[acc2]])

                if args.corr_format == "loadtxt":
                    pdb_mb_dict[pdbid].append((acc1,acc2,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]],is_interaction,surface_area,pdbid,corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)]))
                elif args.corr_format == "read_csv":
                    pdb_mb_dict[pdbid].append((acc2base_species[acc1],acc2base_species[acc2],acc1,acc2,mb_mat[mb_matrix_id_dict[acc1],mb_matrix_id_dict[acc2]],is_interaction,surface_area,pdbid,corr_mat[mb_matrix_uniprot_map_rev[acc1]][mb_matrix_uniprot_map_rev[acc2]]))

    precision, recall, thresholds = precision_recall_curve(interaction_list, correlation_list)
    area = auc(recall, precision)
    print "Pearson PR Area Under Curve: %0.2f" % area
    plt.plot(recall, precision, 'red')

    if nd_matrix_filename != None:
        precision, recall, thresholds = precision_recall_curve(interaction_list, nd_list)
        area = auc(recall, precision)
        print "ND PR Area Under Curve: %0.2f" % area
        plt.plot(recall, precision, 'blue')

    if mb_matrix_filename != None:
        print mb_interaction_list
        print mb_list
        precision, recall, thresholds = precision_recall_curve(mb_interaction_list, mb_list)
        area = auc(recall, precision)
        print "MB stars PR Area Under Curve: %0.2f" % area
        plt.plot(recall, precision, 'black')

    if biogrid_ids_filename != None:
        precision, recall, thresholds = precision_recall_curve(interaction_list, biogrid_list)
        area = auc(recall, precision)
        print "Y2H PR Area Under Curve: %0.2f" % area
        plt.plot(recall, precision, 'green')



    plt.title('Precision-Recall')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    #plt.xlim([0.0, 0.2])
    plt.legend(loc="lower left")

    plt.savefig(args.plot_filename)
    plt.close('all')


    for pdbid in pdb_mb_dict:
        #print pdbid
        for i in pdb_mb_dict[pdbid]:
            #acc1_species, acc2_species, acc1_human,acc2_human,mb_score,is_interaction,surface_area,pdbid,pearson_correlation
            #print "%s %s %s %s %s %s %s %s" % i
            print "%s %s %s %s %s %s %s %s %s" % i


if __name__ == "__main__":
    print "before main"
    main()
