

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

	#kdrew: for every pdb
	for pdbid in pdb_list:

            matrix_ids_filename = args.matrix_ids
            nd_matrix_filename = args.nd_matrix
            corr_matrix_filename = args.corr_matrix

            if args.replace_PDBID:
                matrix_ids_filename = matrix_ids_filename.replace("PDBID",pdbid)
                nd_matrix_filename = nd_matrix_filename.replace("PDBID",pdbid)
                corr_matrix_filename = corr_matrix_filename.replace("PDBID",pdbid)

            #kdrew: read in ids
            matrix_id_file = open(matrix_ids_filename,"rb")
            matrix_id_list = []
            for line in matrix_id_file.readlines():
                    matrix_id_list.append(line.strip())


            nd_mat = np.loadtxt(nd_matrix_filename)
            corr_mat = np.loadtxt(corr_matrix_filename)


            #kdrew: load pisa file
            pisaInt = pdbu.PISA_Interfaces( args.pisa_dir+'/'+pdbid+'_pisa', pdbid=pdbid )

            for acc1, acc2 in it.combinations(matrix_id_list,2):
                surface_area = pisaInt.surface_area_by_acc(acc1, acc2, base_species=args.base_species)
                if surface_area != None:
                        interaction_list.append(1)
                else:
                        interaction_list.append(0)

                nd_list.append( nd_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])
                correlation_list.append( corr_mat[matrix_id_list.index(acc1), matrix_id_list.index(acc2)])

	precision, recall, thresholds = precision_recall_curve(interaction_list, correlation_list)
	area = auc(recall, precision)
	print "PR Area Under Curve: %0.2f" % area
	plt.plot(recall, precision, 'red')

	precision, recall, thresholds = precision_recall_curve(interaction_list, nd_list)
	area = auc(recall, precision)
	print "PR Area Under Curve: %0.2f" % area
	plt.plot(recall, precision, 'blue')

	

	plt.title('Precision-Recall: AUC=%0.2f' % (area,))
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.0])
	plt.legend(loc="lower left")

	plt.savefig(args.plot_filename)
	plt.close('all')


def pdb_calc_interactions(msds, pdbid, base_species, pisa_dir, nd_matrix, nd_proteins):
	print pdbid

	#kdrew: can pass in two pdbids on same line if the complex is split into two pdbs, ex. ribosome 3j3a 3j3b
	split_pdbid = pdbid.split()
	pdbid = split_pdbid[0]
	try:
		pdbid2 = split_pdbid[1]
	except IndexError:
		pdbid2 = None

	#kdrew: for every acc in pdb map to pdb chain
	chain2acc = pu.get_pdb_protein_ids(pdbid)
	if 0 == len(chain2acc):
		print "No accs found for pdbid %s" % (pdbid,)
		return None, None

	#kdrew: for every acc map human version (or specified species)
	acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=base_species)

	#kdrew: get data_set of base_species accs
	data_set = msds.get_data_matrix( )
	acc2idmap = msds.get_id_dict()


	chain_pairs2acc = chain_pairs(msds, chain2acc, acc2base_species)

	a1_dict = calc_correlation( chain_pairs2acc, msds )

	#kdrew: if there is a second pdb in the pdb_list, compile chain combinations within 2nd structure and between 1st and 2nd structure
	if pdbid2 != None:
		chain2acc2 = pu.get_pdb_protein_ids(pdbid2)
		if 0 == len(chain2acc2):
			print "No accs found for pdbid %s" % (pdbid2,)
			return None, None
		acc2base_species2 = pu.get_ortholog( [x for x in chain2acc2.values() if x != None], species1=base_species)
		acc2base_species = dict(acc2base_species.items() + acc2base_species2.items())
		chain_pairs2acc2 = chain_pairs(msds, chain2acc2, acc2base_species)
		a1_dict2 = calc_correlation( chain_pairs2acc2, msds )

		#kdrew: call chain_pairs for inter pdb chain combinations
		chain_pairs2acc_interpdb = chain_pairs(msds, chain2acc, acc2base_species, chain2acc2)
		a1_dict_interpdb = calc_correlation( chain_pairs2acc2, msds )

		#kdrew: combine all correlations into a single dictionary
		a1_dict = dict(a1_dict.items() + a1_dict2.items() + a1_dict_interpdb.items())


        pisaInt = pdbu.PISA_Interfaces( pisa_dir+'/'+pdbid+'_pisa' )
        a2_dict = calc_interactions( chain_pairs2acc, pisaInt )

        nd_mat = np.loadtxt(nd_matrix, delimiter=',')

	#kdrew: convert to data to list format
	a1_list = []
	a2_list = []
	nd_list = []
	for base_key in a1_dict:
		if base_key in a2_dict:
                        try:
                            id1 = nd_proteins.index(base_key[0])
                            id2 = nd_proteins.index(base_key[1])
                        except ValueError:
                            continue
                        print "key: %s corr: %s interaction: %s nd: %s" % (base_key, a1_dict[base_key], a2_dict[base_key], nd_mat[id1,id2])
			a1_list.append( a1_dict[base_key] )
			a2_list.append( a2_dict[base_key] )
                        nd_list.append( nd_mat[id1,id2] ) 


	return (a1_list, a2_list, nd_list)


def chain_pairs(msds, chain2acc, acc2base_species, chain2acc2=None):

	chain_pairs2acc = dict()

	if chain2acc2 == None:
		#kdrew: for every pair of chains 
		#kdrew: find where we have proper mappings to base species acc 
		#kdrew: and acc is in msds data
		for c1, c2 in it.combinations(chain2acc.keys(),2):
			print c1, c2
			base_acc1, base_acc2 = chain_pairs_helper(msds, chain2acc, c1, c2, acc2base_species)

			if base_acc1 == None or base_acc2 == None:
				continue
			else:
				chain_pairs2acc[(c1,c2)] = (base_acc1, base_acc2)
	else:
		for c1 in chain2acc.keys():
			for c2 in chain2acc2.keys():
				print c1, c2
				base_acc1, base_acc2 = chain_pairs_helper(msds, chain2acc, c1, c2, acc2base_species, chain2acc2)

				if base_acc1 == None or base_acc2 == None:
					continue
				else:
					chain_pairs2acc[(c1,c2)] = (base_acc1, base_acc2)


	return chain_pairs2acc

def chain_pairs_helper(msds, chain2acc, c1, c2, acc2base_species, chain2acc2=None):
	acc2idmap = msds.get_id_dict()

	acc1 = chain2acc[c1]
	if chain2acc2 == None:
		acc2 = chain2acc[c2]
	else:
		acc2 = chain2acc2[c2]

	print acc1, acc2

	base_acc1 = None
	base_acc2 = None
	try:
		base_acc1 = acc2base_species[acc1]
	except KeyError:
		print "No mapping to base_species for acc %s" % (acc1)
		return None, None

	try:
		base_acc2 = acc2base_species[acc2]
	except KeyError:
		print "No mapping to base_species for acc %s" % (acc2)
		return None, None

	#kdrew: check to see if accs are in msds
	if base_acc1 not in acc2idmap or base_acc2 not in acc2idmap:
		print "No mapping of base_acc in msds: %s %s" % (base_acc1, base_acc2,)
		return None, None

	return base_acc1, base_acc2

def calc_correlation(chain_pairs2acc, msds):

	data_set = msds.get_data_matrix()
	acc2idmap = msds.get_id_dict()

	corr_dict = dict()
	for c1, c2 in chain_pairs2acc:
		base_acc1 = chain_pairs2acc[(c1,c2)][0]
		base_acc2 = chain_pairs2acc[(c1,c2)][1]

		#kdrew: calculate correlation
		try:
			i = acc2idmap[base_acc1]
		except KeyError:
			print "Missing data in msds: %s" % (base_acc1,)
			continue
		try:
			j = acc2idmap[base_acc2]
		except KeyError:
			print "Missing data in msds: %s" % (base_acc2,)
			continue

		i_set = np.array(data_set[i,])[0]
		j_set = np.array(data_set[j,])[0]
		print "lengths i_set: %s j_set: %s" % (np.count_nonzero(i_set), np.count_nonzero(j_set))
		corr = pearsonr(i_set, j_set)
                print "%s : %s = %s " % (base_acc1, base_acc2, corr[0])
		corr_coeff = corr[0]

		if not np.isnan(corr_coeff): 
			base_key1, base_key2 = base_acc1, base_acc2
			#kdrew: order these based on lex order
			if base_key1 > base_key2:
				base_key1, base_key2 = base_key2, base_key1

			try:
				corr_dict[(base_key1, base_key2)] = min(corr_dict[(base_key1, base_key2)], corr_coeff)
			except KeyError:
				corr_dict[(base_key1, base_key2)] = corr_coeff

	return corr_dict

#kdrew: determine from surface area contact if there is an interaction or not
def calc_interactions(chain_pairs2acc, pisaInt):
	sArea_dict = dict()
	for c1, c2 in chain_pairs2acc:
		base_acc1 = chain_pairs2acc[(c1,c2)][0]
		base_acc2 = chain_pairs2acc[(c1,c2)][1]

		#kdrew: don't calculate between the same protein
		if base_acc1 == base_acc2:
			continue

		sArea = pisaInt.surface_area(c1, c2)

                base_key1, base_key2 = base_acc1, base_acc2
                #kdrew: order these based on lex order
                if base_key1 > base_key2:
                        base_key1, base_key2 = base_key2, base_key1

		if sArea != None:
                        sArea_dict[(base_key1, base_key2)] = 1
                else:
                        sArea_dict[(base_key1, base_key2)] = 0
                        

	return sArea_dict


if __name__ == "__main__":
	main()




