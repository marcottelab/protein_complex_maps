
import numpy as np
import itertools as it
import random 
import argparse
import pickle

import protein_complex_maps.protein_util as pu
import protein_complex_maps.cluster_util as cu
import protein_complex_maps.hierarchical_clustering as hc

def main(): 

	parser = argparse.ArgumentParser(description="Hierarchical clusters fractionation data")
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
	parser.add_argument("--cluster_method", action="store", dest="cluster_method", required=False, default="single",
						help="""Type of linkage clustering, 
						types found: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.cluster.hierarchy.linkage.html, 
						default: single""")

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
	acc2pdbs = pu.map_protein_ids( pdb2acc[args.pdb_id], "ACC", "PDB_ID" )
	#print acc2pdb

	#kdrew: the reverse map has all pdbs with given acc, get only pdb ids with pdb of interest and chain numbering
	acc2pdb = dict()
	pdb_token = args.pdb_id+":"
	for acc in acc2pdbs:
		for id1 in acc2pdbs[acc]:
			if pdb_token in id1:
				try:
					acc2pdb[acc].append(id1)
				except KeyError:
					acc2pdb[acc] = [id1,]

	print acc2pdb


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
					in2pdb[gnames[ortho_mapping[k]]].append(j.split(':')[1])
				except KeyError:
					in2pdb[gnames[ortho_mapping[k]]] = [j.split(':')[1]]

		except KeyError:
			print "No ortholog for %s (%s) : %s" % (k, pdb_gnames[k], acc2pdb[k],)


	interfaces_area = dict()

	#kdrew: get interface area (angstrom^2) from pisa file for every interacting chain
	pisa_file = open(args.pisa_filename, "rb")
	for line in pisa_file.readlines():
		line_split = line.split()
		#interfaces_area[(line_split[1], line_split[6])] = float(line_split[10])
		interfaces_area[(line_split[0], line_split[1])] = float(line_split[2])

	print interfaces_area


	clusters_list = []
	for i in range(len(Y)+1, 2*len(Y)+1):
		clusters_list.append(cu.get_cluster(Y, i))

	gene_clusters_list = []
	for cl in clusters_list:
		gene_cl = []
		for i in cl:
			gene_cl.append(gene_id_map[i])
		gene_clusters_list.append(gene_cl)



	#clusters = []
	#clusters.append(['PSMD6','PSMD3'])
	#clusters.append(['PSMD6','PSMD3','PSMD11'])
	#clusters.append(['PSMD6','PSMD3','PSMD11','PSMD13'])
	#clusters.append(['PSMD6','PSMD3','PSMD11','PSMD13','PSMD12'])
	#clusters.append(['PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7'])
    #
	#clusters.append(['PSMC5','PSMC4'])
    #
	#clusters.append(['PSMD2','PSMC2'])
	#clusters.append(['PSMD2','PSMC2','PSMC1'])
    #
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4'])
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1'])
    #
	#clusters.append(['PSMC6','PSMC3'])
    #
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3'])
    #
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7'])
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7','PSMD14'])
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7','PSMD14','PSMD4'])
    #
    #
	#clusters.append(['PSMA5','PSMA2'])
	#clusters.append(['PSMA5','PSMA2','PSMA6'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3'])
    #
	#clusters.append(['PSMB7','PSMA7'])
    #
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2','PSMA1'])
	#clusters.append(['PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2','PSMA1','PSMA4'])
    #
    #
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7','PSMD14','PSMD4','PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2','PSMA1','PSMA4'])
    #
	#clusters.append(['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7','PSMD14','PSMD4','PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2','PSMA1','PSMA4','PSMD8'])
    #
	#full_complex = ['PSMD2','PSMC2','PSMC1','PSMC5','PSMC4','PSMD1','PSMC6','PSMC3','PSMD6','PSMD3','PSMD11','PSMD13','PSMD12','PSMD7','PSMD14','PSMD4','PSMA5','PSMA2','PSMA6','PSMB3','PSMB1','PSMB4','PSMB6','PSMA3','PSMB7','PSMA7','PSMB5','PSMB2','PSMA1','PSMA4','PSMD8']

	full_complex = gene_id_map.values()

	print clusters_list
	print gene_clusters_list
	#print clusters

	cluster_area = []
	for cl in gene_clusters_list:
		cl_area = calc_cluster_area( cl, in2pdb, interfaces_area )
		cluster_area.append( cl_area )
		#print "%s : %s" % (cl_area, cl,)

		rand_area = []
		for i in range(1000):
			#kdrew: randomly draw len(cl) from full_complex and calc area
			rand_cl_area = calc_cluster_area( random.sample( full_complex, len(cl) ), in2pdb, interfaces_area )
			rand_area.append(rand_cl_area)

		cl_zscore = ( cl_area - np.mean(rand_area) ) / np.std(rand_area)
		print "Iarea: %s Z: %s randmean: %s randstd: %s cluster: %s" % (cl_area, cl_zscore, np.mean(rand_area), np.std(rand_area), cl,)



def calc_cluster_area( cl, pdb_map, interfaces_area ):
	cl_area = 0.0
	#kdrew: for every pair of proteins in cluster
	for i,j in it.combinations(cl,2):
		#print "%s %s %s" % (in2pdb[i], in2pdb[j], interfaces_area[(in2pdb[i],in2pdb[j])], )	
		#kdrew: for every chain of protein1
		for ii in pdb_map[i]:
			#kdrew: for every chain of protein2
			for jj in pdb_map[j]:
				try:
					#kdrew: get interface surface area
					cl_area += interfaces_area[(ii,jj)]
				except KeyError:
					try:
						#print "%s %s %s" % (pdb_map[i], pdb_map[j], interfaces_area[(pdb_map[j],pdb_map[i])], )
						#kdrew: try reverse of chains
						cl_area += interfaces_area[(jj,ii)]
					except KeyError:
						#print "%s %s %s" % (pdb_map[i], pdb_map[j], 0, )
						cl_area += 0.0
	return cl_area



if __name__ == "__main__":
	main()
