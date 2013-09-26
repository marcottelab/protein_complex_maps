
import numpy as np
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc
import pandas as pd

#kdrew: this class is for generating the features of bicluster rows and columns and then creating the logistic regression model

class FeatureGenerator(object):
	row_header = 'row'
	corr_mean_ratio_header = 'corr_mean_ratio'
	intercept_header = 'intercept'
	label_header = 'label'
	
	def __init__(self, bicluster, data_matrix):
		self.__bicluster = bicluster
		self.__data_matrix = data_matrix
		self.__rand_row_bicluster = bc.Bicluster(rows=self.__bicluster.get_random_outside_rows(self.__data_matrix), cols=self.__bicluster.columns() )
		self.__rand_column_bicluster = bc.Bicluster( rows=self.__bicluster.rows(), cols=self.__bicluster.get_random_outside_columns(self.__data_matrix) )

		#kdrew: create a matrix of positives and negatives
		positives = pd.DataFrame( self.__bicluster.rows(), columns=[label_header, row_header, intercept_header,] )
		positives[label_header] = 1
		negatives = pd.DataFrame( self.__rand_row_bicluster.rows(), columns=[label, row_header, intercept_header,] )
		negatives[label_header] = 0
		self.__row_feature_matrix = pd.concat([positives,negatives])



	#kdrew: function to calculate the correlation feature for all rows in and out of bicluster, returns a dataframe of features per row
	def correlation_feature_row(self, ):

		#kdrew: initialize new column and set all values to NAN
		self.__row_feature_matrix[corr_mean_ratio_header] = float('nan')

		#kdrew: first do rows in bicluster and then do rows in random bicluster

		#kdrew: index is the index in the bicluster, i is the index (row number) in the complete data matrix
		for index, i in enumerate(self.__bicluster.rows()):

			print "index:%s, i:%s" % (index, i)
			print self.__bicluster.get_row(self.__data_matrix, row=i)

			#kdrew: get portion of full data matrix that corresponds to the bicluster
			submat = self.__bicluster.get_submatrix( self.__data_matrix )
			print "bicluster submatrix:"
			print submat
			print ""

			#kdrew: calculate correlation of single row vs all other rows
			corrDistIndex = cu.correlation_distribution( submat, index )
			print "correlation_distribution of bicluster submatrix vs index %s:" % (index,)
			print corrDistIndex


			#kdrew: temporarily add row to random bicluster, so we can calculate correlation, remove later
			self.__rand_row_bicluster.add_row(i)

			#kdrew: get portion of full matrix that corresponds to random bicluster
			rand_submat = self.__rand_row_bicluster.get_submatrix( self.__data_matrix )
			print "random bicluster submatrix with row %s:" % (i,)
			print rand_submat


			#kdrew: find index of row number in random bicluster
			#kdrew: ****this makes me a little nervous because we are getting the index from a set which is not ordered****
			rand_bc_index = list(self.__rand_row_bicluster.rows()).index(i)

			#kdrew: calculate correlation of single row vs all other rows in random bicluster
			randbc_corrDistIndex =  cu.correlation_distribution( rand_submat, rand_bc_index )
			print randbc_corrDistIndex

			print corrDistIndex.mean()
			print randbc_corrDistIndex.mean()
			corr_mean_ratio = corrDistIndex.mean()/randbc_corrDistIndex.mean()
			print corr_mean_ratio

			self.__rand_row_bicluster.remove_row(i)

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[row_header]==i, corr_mean_ratio_header ] = corr_mean_ratio




		for index, i in enumerate(self.__rand_row_bicluster.rows()):

			print "index:%s, i:%s" % (index, i)
			print self.__rand_row_bicluster.get_row(self.__data_matrix, row=i)

			#kdrew: get portion of full data matrix that corresponds to the bicluster
			rand_submat = self.__rand_row_bicluster.get_submatrix( self.__data_matrix )
			print "random bicluster submatrix:"
			print rand_submat
			print ""

			#kdrew: calculate correlation of single row vs all other rows
			rand_corrDistIndex = cu.correlation_distribution( rand_submat, index )
			print "correlation_distribution of bicluster submatrix vs index %s:" % (index,)
			print rand_corrDistIndex


			#kdrew: temporarily add row to bicluster, so we can calculate correlation, remove later
			self.__bicluster.add_row(i)

			#kdrew: get portion of full matrix that corresponds to random bicluster
			submat = self.__bicluster.get_submatrix( self.__data_matrix )
			print "bicluster submatrix with row %s:" % (i,)
			print submat


			#kdrew: find index of row number in random bicluster
			#kdrew: ****this makes me a little nervous because we are getting the index from a set which is not ordered****
			bc_index = list(self.__bicluster.rows()).index(i)

			#kdrew: calculate correlation of single row vs all other rows in random bicluster
			bc_corrDistIndex =  cu.correlation_distribution( submat, bc_index )
			print bc_corrDistIndex

			print rand_corrDistIndex.mean()
			print bc_corrDistIndex.mean()
			rand_corr_mean_ratio = rand_corrDistIndex.mean()/bc_corrDistIndex.mean()
			print corr_mean_ratio

			self.__bicluster.remove_row(i)

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[row_header]==i, corr_mean_ratio_header ] = corr_mean_ratio

