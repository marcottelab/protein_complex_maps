
import csv
import pickle
import argparse
import os.path
import sys

import protein_complex_maps.plots.cluster_zscore_vs_dist as czd

def main():

	parser = argparse.ArgumentParser(description="Driver script for linkage clustering of data")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--input_complex_file", action="store", dest="complex_filename", required=False, default=None,
						help="Filename of complex protein identifiers (corum csv)")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--out_pickle_filename", action="store", dest="out_pickle_filename", required=False, default=None,
						help="Filename of output pickle")
	parser.add_argument("--cluster_size_threshold", action="store", type=int, dest="cluster_size_threshold", required=False, default=3,
						help="Size threshold on the number of proteins in a complex, default=3")

	args = parser.parse_args()

	msds = pickle.load( open( args.msds_filename, "rb" ) )

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

	zscore_list_all = []
	distance_list_all = []
	#kdrew: for every complex find some biclusters
	for complex1 in complexes_dict:


		try:
			plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(complex1)+os.path.splitext(args.plot_filename)[1]
			#int_plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(complex1)+'_interactions'+os.path.splitext(args.plot_filename)[1]

			data_set, new_id_map = msds.get_subdata_matrix(complexes_dict[complex1], ignoreNonExistingIds=True)

			#kdrew: do not bother computing small clusters or clusters where we have limited data for
			if new_id_map is None or len(new_id_map) < args.cluster_size_threshold:
				continue

			zscore_list, distance_list = czd.compute_cluster_zscores( complexes_dict[complex1], data_set, new_id_map, msds, cluster_method="single")
			zscore_list_all += zscore_list
			distance_list_all += distance_list
			#plot_cluster_zscore_vs_dist( zscore_list, distance_list, plot_fname)


		#kdrew: if there is a problem continue to next complex
		except:
			e = sys.exc_info()[0]
			print "complex: %s error: %s" % (complex1,e)
			continue

	czd.plot_cluster_zscore_vs_dist( zscore_list_all, distance_list_all, args.plot_filename )

	pickle.dump((zscore_list_all, distance_list_all), open(args.out_pickle_filename,"wb"))


if __name__ == "__main__":
	main()

