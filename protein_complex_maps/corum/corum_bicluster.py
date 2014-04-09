
import csv
import pickle

import protein_complex_maps.coclustering as cc

complex_file = '/home/kdrew/data/corum/allComplexes.csv'
#complex_file = '/home/kdrew/data/corum/fewComplexes.csv'

msds = pickle.load( open( '/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/util/HS_ms2_elutions_msds_peptide_normalized_ids_mapped.p', "rb" ) )

complexes_dict = {}

#kdrew: read in complexes from corum
with open(complex_file, 'rb') as csvfile:
	complexes = csv.reader(csvfile, delimiter=';')
	for row in complexes:
		if row[3] == "Human" and row[4] != '':
			ids = row[4].split(',')
			#kdrew: remove ( ) around extra ids
			complexes_dict[row[0]] = [x.translate(None, '()') for x in ids]


bicluster_dict = {}

#kdrew: for every complex find some biclusters
for complex1 in complexes_dict:
	print "complex id: %s" % (complex1,)

	data_set, new_id_map = msds.get_subdata_matrix(complexes_dict[complex1], ignoreNonExistingIds=True)

	if data_set != None and new_id_map != None:

		#kdrew: run coclustering with different k parameters proportional to 1/2 number of genes
		for k in range(1,int((1.0*len(new_id_map))/2)):
			fit_model = cc.runCluster(data_set, k)
			biclusters = cc.evaluateClusters(msds, data_set, new_id_map, fit_model, None)
			try:
				bicluster_dict[complex1] = bicluster_dict[complex1] + biclusters
			except KeyError:
				bicluster_dict[complex1] = biclusters
				


#kdrew: print out results
for complex1 in bicluster_dict:
	bicluster_dict[complex1].sort(key = lambda tup: tup[1])
	for bicluster in bicluster_dict[complex1]:
		print "complex: %s zscore: %s genes: %s" % (complex1, bicluster[1], ' '.join(bicluster[0]))





