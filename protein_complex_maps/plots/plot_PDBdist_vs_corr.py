

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import argparse
import pickle
import MySQLdb
import itertools as it
from scipy.stats.stats import pearsonr
import Bio.PDB
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
import pandas as pd

import protein_complex_maps.protein_util as pu
import protein_complex_maps.pdb_util as pdbu

def main():

	parser = argparse.ArgumentParser(description="Hierarchical clusters fractionation data")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--pdb_list_filename", action="store", dest="pdb_list_filename", required=True, 
						help="File of benchmark pdbids")
	parser.add_argument("--base_species", action="store", dest="base_species", required=False, default="Hsapiens",
						help="Species to evaluate, default=Hsapiens")
	parser.add_argument("--pickle_results_filename", action="store", dest="pickle_results_filename", required=False, default="./pdb_results_tmp.p",
						help="Filename of pickled results")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--calc_metric", action="store", dest="calc_metric", required=False, default="calc_dist",
						help="""Calculate distance in pdbs ('calc_dist'), requires --pdb_dir 
							or calculate surface area ('calc_surface_area'), requires --pisa_dir
							default 'calc_dist' """)
	parser.add_argument("--pdb_dir", action="store", dest="pdb_dir", required=False, default=None,
						help="Directory of pdbs for calculating distances")
	parser.add_argument("--pisa_dir", action="store", dest="pisa_dir", required=False, default=None,
						help="Directory of pisa files for calculating surface area")

	args = parser.parse_args()

	if args.pdb_dir == None and args.pisa_dir == None:
		print "must set --pdb_dir or --pisa_dir"
		return

	msds = pickle.load( open( args.msds_filename, "rb" ) )

	#kdrew: read in pdb list
	pdb_list_file = open(args.pdb_list_filename,"rb")
	pdb_list = []
	for line in pdb_list_file.readlines():
		pdb_list.append(line.strip())

	pdb_plot_data = dict()

	#kdrew: for every pdb
	for pdbid in pdb_list:
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
			continue


		#kdrew: for every acc map human version (or specified species)
		acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=args.base_species)

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
				continue
			acc2base_species2 = pu.get_ortholog( [x for x in chain2acc2.values() if x != None], species1=args.base_species)
			acc2base_species = dict(acc2base_species.items() + acc2base_species2.items())
			chain_pairs2acc2 = chain_pairs(msds, chain2acc2, acc2base_species)
			a1_dict2 = calc_correlation( chain_pairs2acc2, msds )

			#kdrew: call chain_pairs for inter pdb chain combinations
			chain_pairs2acc_interpdb = chain_pairs(msds, chain2acc, acc2base_species, chain2acc2)
			a1_dict_interpdb = calc_correlation( chain_pairs2acc2, msds )

			#kdrew: combine all correlations into a single dictionary
			a1_dict = dict(a1_dict.items() + a1_dict2.items() + a1_dict_interpdb.items())

		if args.calc_metric == "calc_dist":
			#kdrew: read in pdb
			try:
				structure = Bio.PDB.PDBParser().get_structure(pdbid, args.pdb_dir+'/'+pdbid+'.pdb')
				if pdbid2 != None:
					structure2 = Bio.PDB.PDBParser().get_structure(pdbid2, args.pdb_dir+'/'+pdbid2+'.pdb')
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				continue

			a2_dict = calc_distance( chain_pairs2acc, structure )

			if pdbid2 != None:
				#kdrew: update to deal with interpdb calculations
				a2_dict2 = calc_distance( chain_pairs2acc2, structure2 )
				a2_dict_interpdb = calc_distance( chain_pairs2acc_interpdb, structure, structure2 )
				a2_dict = dict(a2_dict.items() + a2_dict2.items() + a2_dict_interpdb.items())

		elif args.calc_metric == "calc_surface_area":
			pisaInt = pdbu.PISA_Interfaces( args.pisa_dir+'/'+pdbid+'_pisa' )
			a2_dict = calc_surface_area( chain_pairs2acc, pisaInt )

		else:
			print "Error defining calc_metric"
			return


		#kdrew: convert to data to list format
		a1_list = []
		a2_list = []
		for base_key in a1_dict:
			if base_key in a2_dict:
				a1_list.append( a1_dict[base_key] )
				a2_list.append( a2_dict[base_key] )


		pdb_plot_data[pdbid] = (a1_list, a2_list)

	#kdrew: create pickle of current plot data
	pickle.dump(pdb_plot_data, open(args.pickle_results_filename,"wb"))

	plot_data( pdb_plot_data, args.plot_filename )

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
		print corr[0]
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

#kdrew: calculate distance between chains
def calc_surface_area(chain_pairs2acc, pisaInt):
	sArea_dict = dict()
	for c1, c2 in chain_pairs2acc:
		base_acc1 = chain_pairs2acc[(c1,c2)][0]
		base_acc2 = chain_pairs2acc[(c1,c2)][1]

		#kdrew: don't calculate between the same protein
		if base_acc1 == base_acc2:
			continue

		sArea = pisaInt.surface_area(c1, c2)

		if sArea != None:
			base_key1, base_key2 = base_acc1, base_acc2
			#kdrew: order these based on lex order
			if base_key1 > base_key2:
				base_key1, base_key2 = base_key2, base_key1

			#kdrew: if the same acc is represented in multiple chains, take the min dist
			try:
				sArea_dict[(base_key1, base_key2)] = max(sArea_dict[(base_key1, base_key2)], sArea)
			except KeyError:
				sArea_dict[(base_key1, base_key2)] = sArea

	return sArea_dict

#kdrew: calculate distance between chains
def calc_distance(chain_pairs2acc, structure, structure2=None):

	dist_dict = dict()
	for c1, c2 in chain_pairs2acc:
		base_acc1 = chain_pairs2acc[(c1,c2)][0]
		base_acc2 = chain_pairs2acc[(c1,c2)][1]
		#kdrew: don't calculate distance between the same protein
		if base_acc1 == base_acc2:
			continue

		min_dist = pdbu.min_dist(structure, c1, c2, structure2)

		if not np.isnan(min_dist):
			base_key1, base_key2 = base_acc1, base_acc2
			#kdrew: order these based on lex order
			if base_key1 > base_key2:
				base_key1, base_key2 = base_key2, base_key1

			#kdrew: if the same acc is represented in multiple chains, take the min dist
			try:
				dist_dict[(base_key1, base_key2)] = min(dist_dict[(base_key1, base_key2)], min_dist)
			except KeyError:
				dist_dict[(base_key1, base_key2)] = min_dist

	return dist_dict
			
	

def plot_data( pdb_plot_data, plot_filename ):

	total_a1_list = []
	total_a2_list = []

	#kdrew: create a scatter plot with lowess line for each pdbid
	for i, pdbid in enumerate(pdb_plot_data):
		plt.scatter(pdb_plot_data[pdbid][0], pdb_plot_data[pdbid][1])
		ys = lowess(pdb_plot_data[pdbid][1], pdb_plot_data[pdbid][0])
		plt.plot(np.sort(pdb_plot_data[pdbid][0]), ys[:,1], 'red')

		ext = plot_filename.split('.')[-1]
		name = '.'.join(plot_filename.split('.')[:-1])
		pdb_plot_filename =  name + "_" + pdbid + "." + ext

		if plot_filename is None:
			print "plot_filename is None"
			plt.show()
		else:
			plt.savefig(pdb_plot_filename)
			plt.close('all')

		total_a1_list += pdb_plot_data[pdbid][0]
		total_a2_list += pdb_plot_data[pdbid][1]

	#kdrew: scatter plot with lowess line for total set
	plt.scatter(total_a1_list, total_a2_list)
	ys = lowess(total_a2_list, total_a1_list)
	plt.plot(np.sort(total_a1_list), ys[:,1], 'red')
	if plot_filename is None:
		print "plot_filename is None"
		plt.show()
	else:
		plt.savefig(plot_filename)
		plt.close('all')


if __name__ == "__main__":
	main()




