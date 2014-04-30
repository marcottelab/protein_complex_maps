
import csv
import pickle
import argparse

import protein_complex_maps.coclustering as cc
import protein_complex_maps.bicluster.bicluster as bc

def main():

	parser = argparse.ArgumentParser(description="Driver script for coclustering fractionation data")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--iterations", action="store", type=int, dest="iterations", required=False, default=1,
						help="Number of iterations to run, (int)")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, default=None, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--input_complex_file", action="store", dest="complex_filename", required=False, default=None,
						help="Filename of complex protein identifiers (corum csv)")
	parser.add_argument("--scale_threshold", action="store", type=float, dest="scale_threshold", required=False, default=None,
						help="Scale bicluster in data set if less than threshold (float)")
	parser.add_argument("--scale_factor", action="store", type=float, dest="scale_factor", required=False, default=None,
						help="Scale bicluster in data set by given factor (float)")
	#parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
	#					help="Filename of output plot")
	#parser.add_argument("--plot_threshold", action="store", type=float, dest="plot_threshold", required=False, default=0.0,
	#					help="Plot bicluster if score (energy) is less than threshold")

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

	if args.proteins:
		#kdrew: complex does not have an corum id, just mark as input
		complexes_dict['in'] = args.proteins


	bicluster_dict = {}

	#kdrew: for every complex find some biclusters
	for complex1 in complexes_dict:

		data_set, new_id_map = msds.get_subdata_matrix(complexes_dict[complex1], ignoreNonExistingIds=True)
		working_data_set, new_id_map = msds.get_subdata_matrix(complexes_dict[complex1], ignoreNonExistingIds=True)

		for i in xrange(args.iterations):
			print "complex id: %s iteration: %s" % (complex1, i)


			if data_set != None and new_id_map != None:

				#kdrew: run coclustering with different k parameters proportional to 1/2 number of genes
				for k in range(1,int((1.0*len(new_id_map))/2)):
					fit_model = cc.runCluster(working_data_set, k)
					biclusters = cc.evaluateClusters(msds, data_set, new_id_map, fit_model, None)
					try:
						bicluster_dict[complex1] = bicluster_dict[complex1] + biclusters
					except KeyError:
						bicluster_dict[complex1] = biclusters
				

					#kdrew: scale working bicluster
					for bc1 in biclusters:
						if args.scale_threshold and args.scale_factor and bc1[1] < args.scale_threshold:
							bicluster1 = bc.Bicluster(bc1[2], bc1[3])
							working_data_set = bicluster1.scale(working_data_set, args.scale_factor)

	#kdrew: print out results
	for complex1 in bicluster_dict:
		bicluster_dict[complex1].sort(key = lambda tup: tup[1])
		for bicluster in bicluster_dict[complex1]:
			print "complex: %s zscore: %s genes: %s" % (complex1, bicluster[1], ' '.join(bicluster[0]))




if __name__ == "__main__":
	main()

