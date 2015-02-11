

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
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--pdb_list_filename", action="store", dest="pdb_list_filename", required=True, 
						help="File of benchmark pdbids")
	parser.add_argument("--base_species", action="store", dest="base_species", required=False, default="Hsapiens",
						help="Species to evaluate, default=Hsapiens")
	parser.add_argument("--pickle_results_filename", action="store", dest="pickle_results_filename", required=False, default="./interactions_tmp.p",
						help="Filename of pickled results, default = pdb_results_tmp.p")
	parser.add_argument("--load_results_pickle", action="store_true", dest="load_results_pickle", required=False, default=False,
						help="Load results from pickle, overide default filename with --pickle_results_filename")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--pisa_dir", action="store", dest="pisa_dir", required=True, default=None,
						help="Directory of pisa files for calculating surface area")
	parser.add_argument("--nd_proteins", action="store", dest="nd_proteins", nargs='+', required=False, 
						help="Protein ids ordered in ND matrix ")
	parser.add_argument("--nd_matrix", action="store", dest="nd_matrix", required=False, default=None,
						help="Filename of nd matrix")
	args = parser.parse_args()


	#kdrew: read in pdb list
	pdb_list_file = open(args.pdb_list_filename,"rb")
	pdb_list = []
	for line in pdb_list_file.readlines():
		pdb_list.append(line.strip())


	pdb_plot_data = dict()


        a1_list = None
        a2_list = None
	if not args.load_results_pickle:
		msds = pickle.load( open( args.msds_filename, "rb" ) )

		#kdrew: for every pdb
		for pdbid in pdb_list:
			a1_list, a2_list, nd_list = pdb_calc_interactions(msds, pdbid, args.base_species, args.pisa_dir, args.nd_matrix, args.nd_proteins)
			if a1_list != None and a2_list != None:
				pdb_plot_data[pdbid] = (a1_list, a2_list)

		#kdrew: create pickle of current plot data
		pickle.dump(pdb_plot_data, open(args.pickle_results_filename,"wb"))

	else:
		try:
			pdb_plot_data = pickle.load(open(args.pickle_results_filename,"rb"))
		except:
			print "Problem loading results pickle, make sure exists"


        precision, recall, thresholds = precision_recall_curve(a2_list, a1_list)
        area = auc(recall, precision)
        print "PR Area Under Curve: %0.2f" % area
        plt.plot(recall, precision, 'red')

        precision, recall, thresholds = precision_recall_curve(a2_list, nd_list)
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




