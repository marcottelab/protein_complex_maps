
import numpy as np
import itertools as it
import random 
import argparse
import pickle

import protein_complex_maps.protein_util as pu
import protein_complex_maps.cluster_util as cu
import protein_complex_maps.hierarchical_clustering as hc
import protein_complex_maps.pdb_util as pdbu

def main(): 

	parser = argparse.ArgumentParser(description="Hierarchical clusters fractionation data and adds surface area data to each branch")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--pisa_file", action="store", dest="pisa_filename", required=True, 
						help="PDBePISA file containing interface surface area, copy and pasted into file, download options are lacking")
	parser.add_argument("--pdb_id", action="store", dest="pdb_id", required=True, 
						help="PDB id of structure")
	parser.add_argument("--protein_species", action="store", dest="protein_species", required=True, 
						help="species of input proteins")
	parser.add_argument("--pdb_species", action="store", dest="pdb_species", required=True, 
						help="species of pdb id")
	parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
						help="Ignore missing protein ids in msds")
	parser.add_argument("--ortholog_map", action="store", dest="ortholog_map_commandline", nargs='+', required=False, default=[],
						help="Add additional orthologs on commandline when automated ortholog mapping fails (i.e. O15514:P20433 )")
	parser.add_argument("--cluster_method", action="store", dest="cluster_method", required=False, default="single",
						help="""Type of linkage clustering, 
						types found: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.cluster.hierarchy.linkage.html, 
						default: single""")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")

	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )

	data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 

	genename_map = pu.get_genenames_uniprot( new_id_map.values() )
	gene_id_map = dict()
	for i in xrange(len(data_set)):
		#kdrew: sometimes no genename is returned for certain ids, default to original id
		if genename_map[new_id_map[i]] == None:
			gene_id_map[i] = new_id_map[i]
		else:
			gene_id_map[i] = genename_map[new_id_map[i]]

	Y, Y2, D = hc.runCluster( data_set, cluster_method=args.cluster_method )

	#kdrew: got from PDBePISA site and copy and pasted into file, download options are lacking
	#pisa_filename = "4cr2_pisa"
	#pdb_id = "4CR2"

	#kdrew: find all uniprot accs in pdb
	pdb2acc = pu.map_protein_ids( [args.pdb_id,] , "PDB_ID", "ACC" )
	#print pdb2acc
	#kdrew: reverse map uniprot accs to pdb ids
	#acc2pdbs = pu.map_protein_ids( pdb2acc[args.pdb_id], "ACC", "PDB_ID" )
	acc2pdbs = pu.map_protein_ids_to_pdb( pdb2acc[args.pdb_id] )
	#print acc2pdb

	#kdrew: the reverse map has all pdbs with given acc, get only pdb ids with pdb of interest and chain numbering
	acc2pdb = dict()
	#pdb_token = args.pdb_id+":"
	pdb_token = args.pdb_id
	for acc in acc2pdbs:
		for pdbid in acc2pdbs[acc]:
			print "%s %s" % (acc, pdbid,)
			if pdb_token.lower() == pdbid[0].lower():
				pdbid_with_chain = "%s:%s" % (pdbid[0],pdbid[1])
				try:
					acc2pdb[acc].append(pdbid_with_chain)
				except KeyError:
					acc2pdb[acc] = [pdbid_with_chain,]

	print "acc2pdb: %s" % (acc2pdb,)


	#kdrew: map uniprot ids from pdb species to human (or other species of interest)
	#ortho_mapping = pu.get_ortholog( acc2pdb.keys(), 'Hsapiens', 'Scerevisiae' )
	ortho_mapping = pu.get_ortholog( acc2pdb.keys(), args.protein_species, args.pdb_species )
	print ortho_mapping

	print len(acc2pdb)
	print len(ortho_mapping)

	#kdrew: get genenames for all uniprot
	pdb_gnames = pu.get_genenames_uniprot( acc2pdb.keys() )
	gnames = pu.get_genenames_uniprot( ortho_mapping.values() )

	#kdrew: holds a dictionary of lists, key: in_protein -> [pdb chain ids,]
	in2pdb = dict()
	for k in acc2pdb:
		try:
			print " %s (%s) : %s : %s (%s) " % ( k, pdb_gnames[k], acc2pdb[k], ortho_mapping[k], gnames[ortho_mapping[k]], )
			for j in acc2pdb[k]:
				try:
					in2pdb[ortho_mapping[k]].append(j.split(':')[1])
				except KeyError:
					in2pdb[ortho_mapping[k]] = [j.split(':')[1]]

		except KeyError:
			print "No ortholog for %s (%s) : %s, check commandline" % (k, pdb_gnames[k], acc2pdb[k],)
			for j in args.ortholog_map_commandline:
				o1 = j.split(':')[0]
				o2 = j.split(':')[1]
				if k == o1:
					for j in acc2pdb[k]:
						try:
							in2pdb[o2].append(j.split(':')[1])
						except KeyError:
							in2pdb[o2] = [j.split(':')[1]]
				elif k == o2:
					for j in acc2pdb[k]:
						try:
							in2pdb[o1].append(j.split(':')[1])
						except KeyError:
							in2pdb[o1] = [j.split(':')[1]]

	#kdrew: redoing incase of additions from commandline
	gnames = pu.get_genenames_uniprot( in2pdb.keys() )

	#kdrew: get interface area (angstrom^2) from pisa file for every interacting chain
	pisaInt = pdbu.PISA_Interfaces( args.pisa_filename )

	clusters = dict()
	for i in range(len(Y)+1, 2*len(Y)+1):
		clusters[i] = cu.get_cluster(Y, i)

	print "new_id_map: %s" % (new_id_map,)

	gene_clusters = dict()
	for i in clusters:
		gene_cl = []
		for j in clusters[i]:
			gene_cl.append(new_id_map[j])
		gene_clusters[i] = gene_cl


	full_complex = new_id_map.values()

	print clusters
	print gene_clusters

	iarea_gain_list = calc_gene_clusters(gene_clusters, in2pdb, pisaInt, Y, full_complex, gnames)

	rand_iarea_gain_list_total = dict()
	#kdrew: randomize
	for i in xrange(1000):
		#kdrew: randomize new_id_map
		rvals = new_id_map.values()
		random.shuffle(rvals)
		rand_id_map = dict(zip(new_id_map.keys(), rvals))

		rand_gene_clusters = dict()
		for i in clusters:
			gene_cl = []
			for j in clusters[i]:
				gene_cl.append(rand_id_map[j])
			rand_gene_clusters[i] = gene_cl

		rand_iarea_gain_list = calc_gene_clusters(rand_gene_clusters, in2pdb, pisaInt, Y, full_complex, gnames)
		for key in rand_iarea_gain_list:
			try:
				rand_iarea_gain_list_total[key].append(rand_iarea_gain_list[key])
			except KeyError:
				rand_iarea_gain_list_total[key] = [rand_iarea_gain_list[key]]


	rand_iarea_gain_list_mean = dict()
	rand_iarea_gain_list_std = dict()
	for key in rand_iarea_gain_list_total:
		rand_iarea_gain_list_mean[key] = np.mean(rand_iarea_gain_list_total[key])
		rand_iarea_gain_list_std[key] = np.std(rand_iarea_gain_list_total[key])


	if args.plot_filename:
 		circle_annotate = [iarea_gain_list[cl] for cl in sorted(iarea_gain_list.keys()) ] 
 		circle_annotate_rand = [rand_iarea_gain_list_mean[cl] for cl in sorted(rand_iarea_gain_list_mean.keys()) ] 
		print circle_annotate
		print "monotonic? %s" % (Y[:,2], )
		hc.plotJustDendrogram(Y, args.plot_filename, gene_id_map, circle_annotate =  circle_annotate )
		hc.plotJustDendrogram(Y, args.plot_filename+".rand.pdf", gene_id_map, circle_annotate =  circle_annotate_rand )
	

def calc_gene_clusters(gene_clusters, in2pdb, pisaInt, Y, full_complex, gnames):
	iarea_gain_list = dict()
	for cl in gene_clusters:
		cl_area = calc_cluster_area( gene_clusters[cl], in2pdb, pisaInt )
		ch1, ch2 = cu.get_cluster_children(Y, cl)
		#print "ch1: %s ch2: %s" % (ch1, ch2,)
		if ch1 > len(Y):
			ch1_area = calc_cluster_area( gene_clusters[ch1], in2pdb, pisaInt )
		else:
			ch1_area = 0.0

		if ch2 > len(Y):
			ch2_area = calc_cluster_area( gene_clusters[ch2], in2pdb, pisaInt )
		else:
			ch2_area = 0.0

		#print "%s : %s" % (cl_area, cl,)

		#rand_area = []
		#for i in range(1000):
		#	#kdrew: randomly draw len(gene_clusters[cl]) from full_complex and calc area
		#	rand_cl_area = calc_cluster_area( random.sample( full_complex, len(gene_clusters[cl]) ), in2pdb, pisaInt )
		#	rand_area.append(rand_cl_area)
        #
		#cl_zscore = ( cl_area - np.mean(rand_area) ) / np.std(rand_area)
		#print "Total Iarea: %s Iarea_gain: %s Z: %s randmean: %s randstd: %s cluster: %s" % (cl_area, cl_area - ch1_area - ch2_area, cl_zscore, np.mean(rand_area), np.std(rand_area), [gnames[i] for i in gene_clusters[cl]],)
		print "Total Iarea: %s Iarea_gain: %s cluster: %s" % (cl_area, cl_area - ch1_area - ch2_area, [gnames[i] for i in gene_clusters[cl]],)
		iarea_gain_list[cl] = cl_area - ch1_area - ch2_area
	
	return iarea_gain_list




def calc_cluster_area( cl, pdb_map, pisaInt ):
	cl_area = 0.0
	#kdrew: for every pair of proteins in cluster
	for i,j in it.combinations(cl,2):
		#print "%s %s %s" % (in2pdb[i], in2pdb[j], pisaInt[(in2pdb[i],in2pdb[j])], )	
		#kdrew: for every chain of protein1
		for ii in pdb_map[i]:
			#kdrew: for every chain of protein2
			for jj in pdb_map[j]:
				#kdrew: get interface surface area
				sa = pisaInt.surface_area(ii,jj)
				if None != sa:
					cl_area += sa
				else:
					cl_area += 0.0
	return cl_area



if __name__ == "__main__":
	main()
