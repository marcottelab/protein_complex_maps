
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
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
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


	msds = pickle.load( open( args.msds_filename, "rb" ) )

	data_set, new_id_map = msds.get_subdata_matrix(protein_ids, ignoreNonExistingIds=True)

	out_id_file = open(os.path.splitext(args.out_filename)[0]+'.ids', "wb")
	for i in xrange(len(new_id_map)):
		out_id_file.write(new_id_map[i]+'\n')

	out_id_file.close()
	matrix_clean = np.nan_to_num(data_set)
	corrcoefMat = np.nan_to_num(np.corrcoef(matrix_clean))

	np.savetxt(args.out_filename, corrcoefMat)



if __name__ == "__main__":
	main()


