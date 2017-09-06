

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
import Bio.PDB
import matplotlib.pyplot as plt

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

import protein_complex_maps.protein_util as pu
import protein_complex_maps.pdb_util as pdbu

def main():

	parser = argparse.ArgumentParser(description="Plot precision recall of protein interactions for network deconvolution and correlation")
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
	parser.add_argument("--corr_matrix", action="store", dest="corr_matrix", required=False, default=None,
						help="Filename of corr matrix")
	parser.add_argument("--matrix_ids", action="store", dest="matrix_ids", required=False, default=None,
						help="Filename of ids in matrix")
	parser.add_argument("--biogrid_matrix", action="store", dest="biogrid_matrix", required=False, default=None,
						help="Filename of biogrid matrix")
	parser.add_argument("--biogrid_ids", action="store", dest="biogrid_ids", required=False, default=None,
						help="Filename of ids in matrix from biogrid")
	parser.add_argument("--replace_PDBID", action="store_true", dest="replace_PDBID", required=False, default=False,
						help="Flag to replace PDBID in nd_matrix, corr_matrix and matrix_ids args with pdb ids in pdb_list_filename")
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
	correlation_list = []
	nd_list = []
        biogrid_list = []

	#kdrew: for every pdb
	for pdbid in pdb_list:

            chain2acc = pu.get_pdb_protein_ids(pdbid)
            acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=args.base_species, reversible=True)

            matrix_ids_filename = args.matrix_ids
            biogrid_ids_filename = args.biogrid_ids
            nd_matrix_filename = args.nd_matrix
            corr_matrix_filename = args.corr_matrix
            biogrid_matrix_filename = args.biogrid_matrix

            if args.replace_PDBID:
                matrix_ids_filename = matrix_ids_filename.replace("PDBID",pdbid)
                biogrid_ids_filename = biogrid_ids_filename.replace("PDBID",pdbid)
                nd_matrix_filename = nd_matrix_filename.replace("PDBID",pdbid)
                corr_matrix_filename = corr_matrix_filename.replace("PDBID",pdbid)
                biogrid_matrix_filename = biogrid_matrix_filename.replace("PDBID",pdbid)

            #kdrew: read in ids
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

            #kdrew: readin matrices
            if nd_matrix_filename != None:
                nd_mat = np.loadtxt(nd_matrix_filename)
            corr_mat = np.loadtxt(corr_matrix_filename)

            if biogrid_ids_filename != None:
                biogrid_mat = np.loadtxt(biogrid_matrix_filename)


            #kdrew: load pisa file
            pisaInt = pdbu.PISA_Interfaces( args.pisa_dir+'/'+pdbid+'_pisa', pdbid=pdbid )

            for acc1, acc2 in it.combinations(pdb_matrix_id_list,2):
                surface_area = pisaInt.surface_area_by_acc(acc1, acc2, base_species=args.base_species)
                if surface_area != None:
                        interaction_list.append(1)
                else:
                        interaction_list.append(0)

                if nd_matrix_filename != None:
                    nd_list.append( nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])

                correlation_list.append( corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])

                if biogrid_ids_filename != None:
                    biogrid_list.append( biogrid_mat[biogrid_id_list.index(acc1), biogrid_id_list.index(acc2)])

	precision, recall, thresholds = precision_recall_curve(interaction_list, correlation_list)
	area = auc(recall, precision)
	print "PR Area Under Curve: %0.2f" % area
	plt.plot(recall, precision, 'red')

        if nd_matrix_filename != None:
            precision, recall, thresholds = precision_recall_curve(interaction_list, nd_list)
            area = auc(recall, precision)
            print "PR Area Under Curve: %0.2f" % area
            plt.plot(recall, precision, 'blue')

        if biogrid_ids_filename != None:
            precision, recall, thresholds = precision_recall_curve(interaction_list, biogrid_list)
            area = auc(recall, precision)
            print "PR Area Under Curve: %0.2f" % area
            plt.plot(recall, precision, 'green')

	

	#plt.title('Precision-Recall: AUC=%0.2f' % (area,))
	plt.title('Precision-Recall')
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	#plt.xlim([0.0, 0.2])
	plt.legend(loc="lower left")

	plt.savefig(args.plot_filename)
	plt.close('all')




if __name__ == "__main__":
	main()




