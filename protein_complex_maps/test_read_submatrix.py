
import numpy as np
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc

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

bicluster1 = bc.Bicluster(rows = [3,10,13,14,22,23,26], cols = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70])

bc1_submat = bicluster1.get_submatrix(clean_data_matrix)
print bc1_submat

#bc1_corrDist = cu.correlation_distribution(bc1_submat)
#print bc1_corrDist



