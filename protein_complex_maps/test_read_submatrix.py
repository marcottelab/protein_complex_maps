
import numpy as np
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc
from scipy.stats import ks_2samp

sample_filename = "/Users/kdrew/data/bborgeson/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"

sample_file = open(sample_filename, 'rb')
#kdrew: eat header
line = sample_file.readline()

data = []

for line in sample_file.readlines():
	#print line
	line_data = line.split()
	line_array = map(float,line_data[2:])
	print line_array
	data.append(line_array)

#print data

data_matrix = np.asmatrix(data)

#print data_matrix


#kdrew: remove columns with all zeros
clean_data_matrix = data_matrix.compress(~np.array(np.all(data_matrix[:]==0,axis=0))[0],axis=1)
print clean_data_matrix
##kdrew: remove columns with all zeros
#self.__clean_data_matrix = self.__data_matrix.compress(~np.array(np.all(self.__data_matrix[:]==0,axis=0))[0],axis=1)

bicluster1 = bc.Bicluster(rows = [3,10,13,14,22,26], cols = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70])

bc1_submat = bicluster1.get_submatrix(clean_data_matrix)
print bc1_submat

bc1_corrDist = cu.correlation_distribution(bc1_submat)
print bc1_corrDist

bc1_r0_corrDist = cu.correlation_distribution(bc1_submat, 0)
print bc1_r0_corrDist

rand_bicluster1 = bc.Bicluster(rows=bicluster1.get_random_outside_rows(clean_data_matrix), cols=bicluster1.columns())

#kdrew: index is the index in the bicluster, i is the index (row number) in the complete data matrix
for index, i in enumerate(bicluster1.rows()):
	print "index:%s, i:%s" % (index, i)
	print bicluster1.get_row(clean_data_matrix, row=i)
	submat = bicluster1.get_submatrix( clean_data_matrix )
	print "bicluster1 submatrix:"
	print submat

	print ""
	submatWO = bicluster1.get_submatrix( clean_data_matrix, without_rows=[i] )
	corrDistWO = cu.correlation_distribution( submatWO )
	print "correlation_distribution without %s:" % (i,)
	print corrDistWO

	corrDistIndex = cu.correlation_distribution( submat, index )
	print "correlation_distribution of bicluster submatrix vs index %s:" % (index,)
	print corrDistIndex

	#ks_biclust = ks_2samp(corrDistWO, corrDistIndex)
	#print ks_biclust

	rand_bicluster1.add_row(i)

	rand_submat = rand_bicluster1.get_submatrix( clean_data_matrix )
	print "random bicluster submatrix with row %s:" % (i,)
	print rand_submat
	rand_submatWO = rand_bicluster1.get_submatrix( clean_data_matrix, without_rows=[i])

	randbc_corrDistWO = cu.correlation_distribution( rand_submatWO )
	print randbc_corrDistWO

	rand_bc_index = list(rand_bicluster1.rows()).index(i)
	randbc_corrDistIndex =  cu.correlation_distribution( rand_submat, rand_bc_index )
	print randbc_corrDistIndex

	#ks_rand_biclust = ks_2samp(randbc_corrDistWO, randbc_corrDistIndex)
	#print ks_rand_biclust

	#kdrew: compute ratio of ks tests, this gives a value of how correlated the bicluster is compared to the background
	#kdrew: values of ~1.0 suggest no correlation above background, small values suggest bicluster is more correlated than background
	#print ks_biclust[0]/ks_rand_biclust[0]

	print corrDistIndex.mean()
	print randbc_corrDistIndex.mean()
	print corrDistIndex.mean()/randbc_corrDistIndex.mean()

	print ""

	rand_bicluster1.remove_row(i)




