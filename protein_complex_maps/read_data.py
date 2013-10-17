
import numpy as np


def read_datafile(fhandle, header=True):
	if header:
		#kdrew: eat header
		line = fhandle.readline()
	
	data = []
	name_list = []

	for line in fhandle.readlines():
		#print line
		line_data = line.split()
		name_list.append(line_data[0])
		line_array = map(float,line_data[2:])
		print line_array
		data.append(line_array)

	#print data

	data_matrix = np.asmatrix(data)

	return data_matrix, name_list

