
import logging
import numpy as np

import protein_complex_maps.normalization_util as nu
import protein_complex_maps.protein_util as pu

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s')

class MSDataSet(object):

	def __init__( self ):
		self.__master_data_matrix = None
		#kdrew: holds original ids read from file and used to link other files read in, ordered by matrix rows
		self.__master_name_list = None
		self.__master_fraction_list = None
		#kdrew: holds mapping of protein ids to matrix indices, includes ids from master_name_list
		self.__id_dict = dict()
		self.__frac_dict = dict()

	def load_file( self, file_handle, header=False, normalize=False):
		
		data_matrix1, name_list1 = read_datafile(file_handle, header=header)
		if normalize:
			data_matrix1 = nu.normalize_over_columns(data_matrix1)
		if self.__master_data_matrix == None:
			self.__master_data_matrix = data_matrix1
			self.__master_name_list = name_list1
		else:
			self.__master_data_matrix, self.__master_name_list = concat_data_matrix( self.__master_data_matrix, self.__master_name_list, data_matrix1, name_list1)

		self.update_id_dict()

	def update_id_dict(self, reset=False):
		if reset:
			self.__id_dict = dict()
			self.__frac_dict = dict()
		#kdrew: updating id_dict with current name list
		for i, name in enumerate(self.__master_name_list):
			self.__id_dict[name] = i

		if self.__master_fraction_list != None:
			for i, frac in enumerate(self.__master_fraction_list):
				self.__frac_dict[frac] = i

	#kdrew: populate msds with peptide count dictionary
	#kdrew: dictionary should be of the form protein->fraction->peptide->count
	def create_by_peptide_counts( self, protein_counts ):
		#kdrew: get a list of all proteins
		protein_list = protein_counts.keys()

		#kdrew: get a list of all fractions
		fractions = set()
		for prot in protein_list:
			fractions = fractions.union(protein_counts[prot].keys())
		fractions_list = list(fractions)
		fractions_list.sort()
		#print fractions_list

		#kdrew: create data matrix size: proteins X fractions
		dmat = np.matrix(np.zeros(shape=(len(protein_list),len(fractions_list))))

		for i, prot in enumerate(protein_list):
			for j, fraction in enumerate(fractions_list):
				try:
					 for peptide in protein_counts[prot][fraction]:
						#kdrew: here is where normalization should be done
						dmat[i,j] += protein_counts[prot][fraction][peptide]
				except KeyError:
					continue
	
		#print dmat

		if self.__master_data_matrix == None:
			self.__master_data_matrix = dmat
			self.__master_name_list = protein_list
			self.__master_fraction_list = fractions_list
		else:
			self.__master_data_matrix, self.__master_name_list = concat_data_matrix( self.__master_data_matrix, self.__master_name_list, data_matrix1, name_list1)
			self.__master_fraction_list += fractions_list

		self.update_id_dict()


	#kdrew: adds mappings of protein ids to matrix indices
	def map_ids( self, from_id, to_id):
		#kdrew: map master_name_list from current db_id to db_id
		#kdrew: update master_name_list and current db_id
		protids_map = pu.map_protein_ids( self.__master_name_list, from_id, to_id )
		for protid in protids_map:
			i = self.__master_name_list.index(protid)
			for mapped_id in protids_map[protid]:
				self.__id_dict[mapped_id] = i
		
		#return self.__id_dict

	#kdrew: dictionary of protein ids (keys) to matrix indices
	def get_id_dict( self ):
		return self.__id_dict

	def get_fraction_dict( self ):
		return self.__frac_dict


	def set_id_dict( self, id_dict ):
		self.__id_dict = id_dict

	#def get_data_matrix( self, names=None, remove_zero=False ):

	#kdrew: get submatrix based on protein identifiers
	def get_subdata_matrix( self, ids, ignoreNonExistingIds=False):
		matrix = self.get_data_matrix()
		id_indices = []
		new_map = dict()
		for i,i_d in enumerate(ids):
			#print i, i_d
			try:
				index1 = self.__id_dict[i_d]
				if index1 not in id_indices:
					id_indices.append(index1)
					new_map[i] = i_d
			except KeyError:
				if ignoreNonExistingIds:
					print "get_subdata_matrix ignoring: %s" % (i_d,)
					continue
				else:
					raise


		#print "id_indices: %s mat_shape: %s" % (id_indices, matrix.shape[1], )
		assert matrix.shape[1] > 0, "no columns in matrix" 
		if len(id_indices) > 0:
			#kdrew: set ensures no duplicates, but I lose mapping
			return matrix[np.ix_(id_indices, range(0,matrix.shape[1]))], new_map
		else:
			return None, None
		
	def get_name2index( self, ):
		return {v:k for k, v in self.__id_dict.items()}


	def get_data_matrix( self, remove_zero=False ):
		#print "get_data_matrix"
		#print self.__master_data_matrix

		##kdrew: might want to return new msds object because can no longer map protein ids to indices of new matrix
		#if names != None:
		#	rows = []
		#	for name in names:
		#		#print "name: %s, index: %s" % (name, self.__master_name_list.index(name))
		#		rows.append(self.__id_dict[name]))
        #
		#	cols = range(0,self.__master_data_matrix.shape[1])
		#	submatrix = self.__master_data_matrix[np.ix_(rows, cols)]
        #
		#	#kdrew: only remove zero columns
		#	if remove_zero:
		#		submatrix = nu.remove_zero(submatrix, zero_rows=False, zero_columns=True)
        #
		#	return submatrix
        
		if remove_zero:
			return nu.remove_zero(self.__master_data_matrix, zero_rows=False, zero_columns=True)

		return self.__master_data_matrix

	def set_data_matrix( self, data_matrix ):
		self.__master_data_matrix = data_matrix


def read_datafile(fhandle, header=True):
	if header:
		#kdrew: eat header
		line = fhandle.readline()
	
	data = []
	name_list = []

	for line in fhandle.readlines():
		print line
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
	



