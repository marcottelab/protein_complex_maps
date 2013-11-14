import logging
import numpy as np

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s')

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
		logging.debug(line_array)
		data.append(line_array)

	#print data

	data_matrix = np.asmatrix(data)

	return data_matrix, name_list

def concat_data_matrix(data_matrix1, name_list1, data_matrix2, name_list2):
	
	name_set1 = set(name_list1)
	name_set2 = set(name_list2)

	dm1_num_columns = data_matrix1.shape[1]
	dm2_num_columns = data_matrix2.shape[1]

	combined_set = name_set1.union(name_set2)
	logging.debug(combined_set)

	return_mat = None
	return_name_list = []
                        
	for name in combined_set:
		logging.debug("name: %s" % (name,))
		return_name_list.append(name)
		try:
			idx1 = name_list1.index(name)
			row1 = data_matrix1[idx1]
			logging.debug("idx1: %s" % (idx1,))
			logging.debug("row1: %s" % (row1,))
		except ValueError:
			#kdrew: not in list
			row1 = np.zeros(dm1_num_columns)

		try:
			idx2 = name_list2.index(name)
			row2 = data_matrix2[idx2]
		except ValueError:
			#kdrew: not in list
			row2 = np.zeros(dm2_num_columns)

		complete_row = np.append(np.array(row1), np.array(row2))
	

		#kdrew: add complete row to full matrix
		if return_mat == None:
			return_mat = np.matrix(complete_row)
		else:
			return_mat = np.vstack([return_mat,complete_row])

	return return_mat, return_name_list
	



