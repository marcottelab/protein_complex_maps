
import csv
import pickle
import argparse
import os.path
import sys


import protein_complex_maps.plots.cluster_zscore_vs_dist as czd
import protein_complex_maps.physical_interaction_clustering as pic
import protein_complex_maps.protein_util as pu
import itertools as it

def main():

	parser = argparse.ArgumentParser(description="Driver script for linkage clustering of data")
	#parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
	#					help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--input_complex_file", action="store", dest="complex_filename", required=False, default=None,
						help="Filename of complex protein identifiers (corum csv)")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--out_pickle_filename", action="store", dest="out_pickle_filename", required=False, default=None,
						help="Filename of output pickle")
	parser.add_argument("--cluster_size_threshold", action="store", type=int, dest="cluster_size_threshold", required=False, default=3,
						help="Size threshold on the number of proteins in a complex, default=3")

	args = parser.parse_args()

	#msds = pickle.load( open( args.msds_filename, "rb" ) )

	complexes_dict = {}

	if args.complex_filename: 
		#kdrew: read in complexes from corum
		with open(args.complex_filename, 'rb') as csvfile:
			complexes = csv.reader(csvfile, delimiter=';')
			for row in complexes:
				if row[3] == "Human" and row[4] != '':
					ids = row[4].split(',')
					#kdrew: remove ( ) around extra ids
					complexes_dict[row[0]] = [x.translate(None, '()') for x in ids]

	#kdrew: convert ids to biogrid and search biogrid for interactions
	biogrid_map = pu.map_protein_ids( list(it.chain(*complexes_dict.values()))  , "ACC+ID", 'BIOGRID_ID' )

	#kdrew: for every complex find some biclusters
	for complex1 in complexes_dict:

		try:
			#plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(complex1)+os.path.splitext(args.plot_filename)[1]
			#int_plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(complex1)+'_interactions'+os.path.splitext(args.plot_filename)[1]

			#data_set, new_id_map = msds.get_subdata_matrix(complexes_dict[complex1], ignoreNonExistingIds=True)

			#kdrew: do not bother computing small clusters or clusters where we have limited data for
			if len(complexes_dict[complex1]) < args.cluster_size_threshold:
				continue

			D = pic.create_interaction_matrix( complexes_dict[complex1], biogrid_map, experimental_system = ['Co-crystal Structure','Reconstituted Complex'] )
			print "complex: %s proteins: %s interactions: %s" % ( complex1, len(complexes_dict[complex1]), D.sum() )

		#kdrew: if there is a problem continue to next complex
		except:
			e = sys.exc_info()[0]
			print "complex: %s error: %s" % (complex1,e)
			continue


if __name__ == "__main__":
	main()

