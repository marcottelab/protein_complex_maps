
import os.path 

import numpy as np

import argparse
import pickle
import itertools as it
import Bio.PDB

import protein_complex_maps.protein_util as pu
import protein_complex_maps.pdb_util as pdbu

def main():

	parser = argparse.ArgumentParser(description="Output correlation matrix")
	parser.add_argument("--biogrid_filename", action="store", dest="biogrid_filename", required=True, 
						help="Filename of BioGrid interactions")
	parser.add_argument("--pdb_list_filename", action="store", dest="pdb_list_filename", required=False, default=None,
						help="File of benchmark pdbids")
	parser.add_argument("--pdbid", action="store", dest="pdbid", required=False, default=None, 
						help="PDB id")
	parser.add_argument("--base_species", action="store", dest="base_species", required=False, default="Hsapiens",
						help="Species to evaluate, default=Hsapiens")
	parser.add_argument("--out_filename", action="store", dest="out_filename", required=False, default="matrix.txt",
						help="Name of output matrix, ids will be stored in 'out_filename'.ids")
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

	protein_ids = []
	for pdbid in pdb_list:
		chain2acc = pu.get_pdb_protein_ids(pdbid)
		acc2base_species = pu.get_ortholog( [x for x in chain2acc.values() if x != None], species1=args.base_species)
		protein_ids = protein_ids + acc2base_species.values()


        #kdrew: map protein_ids to biogrid ids
        biogrid_map = pu.map_protein_ids(protein_ids, 'ACC', 'BIOGRID_ID')
        print biogrid_map

        #kdrew: collapse list to single entry
        biogrid_map = {k:biogrid_map[k][0] for k in biogrid_map if len(biogrid_map[k]) > 0}
        print biogrid_map

        #kdrew: make map reversible
        biogrid_map.update({biogrid_map[x]:x for x in biogrid_map})

        biogrid_interaction_dict = dict()
        #kdrew: readin biogrid file
        bg_file = open(args.biogrid_filename,"rb")
        for line in bg_file.readlines():
            if line.split()[0] != "#BioGRID":
                biogrid_idA = line.split()[3]
                biogrid_idB = line.split()[4]
                #print biogrid_idA
                #print biogrid_idB
                if biogrid_idA in biogrid_map.values() and biogrid_idB in biogrid_map.values():
                    biogrid_interaction_dict[biogrid_idA] = biogrid_idB
                    biogrid_interaction_dict[biogrid_idB] = biogrid_idA
                    print line
    
        print biogrid_interaction_dict

        data_mat = np.zeros([len(protein_ids),len(protein_ids)])
        for k in biogrid_interaction_dict:
            #kdrew: put a 1.0 for every interaction
            data_mat[protein_ids.index(biogrid_map[k]),protein_ids.index(biogrid_map[biogrid_interaction_dict[k]])] = 1.0

        for k in protein_ids:
            #kdrew: put a 1.0 for every self
            data_mat[protein_ids.index(k),protein_ids.index(k)] = 1.0

        print data_mat

	out_id_file = open(os.path.splitext(args.out_filename)[0]+'.ids', "wb")
        for i in protein_ids:
		out_id_file.write(i+'\n')
	out_id_file.close()

	np.savetxt(args.out_filename, data_mat)



if __name__ == "__main__":
	main()


