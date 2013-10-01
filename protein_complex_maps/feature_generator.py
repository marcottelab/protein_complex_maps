
import numpy as np
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc
import pandas as pd

#kdrew: this class is for generating the features of bicluster rows and columns and then creating the logistic regression model

class FeatureGenerator(object):
	#kdrew: define labels in pandas dataframe
	column_header = 'column'
	row_header = 'row'
	corr_mean_ratio_header = 'corr_mean_ratio'
	corr_mean_bc_header = 'corr_mean_bc'
	corr_mean_rand_header = 'corr_mean_rand'
	corr_gain_bc_header = 'corr_gain_bc'
	corr_gain_rand_header = 'corr_gain_rand'
	intercept_header = 'intercept'
	label_header = 'label'
	
	def __init__(self, bicluster, data_matrix, seed=None):
		self.__bicluster = bicluster
		self.__data_matrix = data_matrix
		self.__rand_row_bicluster = bc.Bicluster(rows=self.__bicluster.get_random_outside_rows(self.__data_matrix,seed=seed), cols=self.__bicluster.columns() )
		self.__rand_column_bicluster = bc.Bicluster( rows=self.__bicluster.rows(), cols=self.__bicluster.get_random_outside_cols(self.__data_matrix, seed=seed ) )

		#kdrew: create a dataframe of positive and negative items to feed into regression
		row_positives = pd.DataFrame( self.__bicluster.rows(), columns=[self.row_header]) 
		row_positives[self.label_header] = 1
		row_negatives = pd.DataFrame( self.__rand_row_bicluster.rows(), columns=[self.row_header]) 
		row_negatives[self.label_header] = 0
		self.__row_feature_matrix = pd.concat([row_positives,row_negatives], ignore_index=True)

		#kdrew: this might get a little confusing because we are now building a feature matrix on bicluster columns 
		#kdrew: which will be rows in the feature matrix
		column_positives = pd.DataFrame( self.__bicluster.columns(), columns=[self.column_header] )
		column_positives[self.label_header] = 1
		column_negatives = pd.DataFrame( self.__rand_column_bicluster.columns(), columns=[self.column_header] )
		column_negatives[self.label_header] = 0
		self.__column_feature_matrix = pd.concat([column_positives,column_negatives], ignore_index=True)


	def get_row_feature_matrix(self,):
		return self.__row_feature_matrix

	def get_column_feature_matrix(self, ):
		return self.__column_feature_matrix

	#kdrew: function to calculate the correlation feature for all columns in and out of bicluster, 
	#kdrew: creates a dataframe of features per bicluster column
	def correlation_feature_column(self, ):

		#kdrew: initialize new column and set all values to NAN
		self.__column_feature_matrix[self.corr_gain_bc_header] = float('nan')
		self.__column_feature_matrix[self.corr_gain_rand_header] = float('nan')

		print self.__column_feature_matrix

		#kdrew: first do columns in bicluster and then do columns in random bicluster

		#kdrew: index is the index in the bicluster, i is the index (column number) in the complete data matrix
		for index, i in enumerate(self.__bicluster.columns()):

			print "index: %s, i: %s" % (index, i)
			corr_gain, rand_corr_gain = self.correlation_feature_column_by_index( self.__bicluster, self.__rand_column_bicluster, index, i )

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where column_header == i and column is corr_mean_ratio_header
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_bc_header ] = corr_gain
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_rand_header ] = rand_corr_gain

		for index, i in enumerate(self.__rand_column_bicluster.columns()):

			rand_corr_gain, bc_corr_gain = self.correlation_feature_column_by_index( self.__rand_column_bicluster, self.__bicluster, index, i )

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all columns where column_header == i and column is corr_mean_ratio_header
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_bc_header ] = bc_corr_gain
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_rand_header ] = rand_corr_gain

		print self.__column_feature_matrix


	#kdrew: function to calculate the correlation feature for all rows in and out of bicluster, 
	#kdrew: creates a dataframe of features per row
	def correlation_feature_row(self, ):

		#kdrew: initialize new column and set all values to NAN
		self.__row_feature_matrix[self.corr_mean_bc_header] = float('nan')
		self.__row_feature_matrix[self.corr_mean_rand_header] = float('nan')
		self.__row_feature_matrix[self.corr_mean_ratio_header] = float('nan')

		#kdrew: first do rows in bicluster and then do rows in random bicluster

		#kdrew: index is the index in the bicluster, i is the index (row number) in the complete data matrix
		for index, i in enumerate(self.__bicluster.rows()):

			corr_mean, rand_corr_mean = self.correlation_feature_row_by_index( self.__bicluster, self.__rand_row_bicluster, index, i )
			#print corrDistIndex.mean()
			#print bc2_corrDistIndex.mean()
			corr_mean_ratio = corr_mean/rand_corr_mean
			#print corr_mean_ratio


			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_bc_header ] = corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_rand_header ] = rand_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_ratio_header ] = corr_mean_ratio

		for index, i in enumerate(self.__rand_row_bicluster.rows()):

			rand_corr_mean, bc_corr_mean = self.correlation_feature_row_by_index( self.__rand_row_bicluster, self.__bicluster, index, i )

			rand_corr_mean_ratio = bc_corr_mean/rand_corr_mean

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_bc_header ] = bc_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_rand_header ] = rand_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_ratio_header ] = rand_corr_mean_ratio

		#print self.__row_feature_matrix

	#kdrew: when calling this on columns in bicluster, bicluster1 = bicluster, bicluster2 = random bicluster
	#kdrew: switch them when calling on columns from random bicluster
	def correlation_feature_column_by_index( self, bicluster1, bicluster2, index, i ):

			#kdrew: get portion of full data matrix that corresponds to the bicluster
			submat = bicluster1.get_submatrix( self.__data_matrix )
			#kdrew: calculate correlation of all rows in bicluster
			corrDist = cu.correlation_distribution( submat )

			#kdrew: remove column i from bicluster, regenerate submat
			bicluster1.remove_column( i )
			submat = bicluster1.get_submatrix( self.__data_matrix )
			#kdrew: recalculate correlation of all rows in bicluster without column i
			corrDist_without_i = cu.correlation_distribution( submat )

			#kdrew: re-add column i to bicluster
			bicluster1.add_column( i )
			bc2_submat = bicluster2.get_submatrix( self.__data_matrix )
			#kdrew: calculate correlation of all rows in bicluster2
			bc2_corrDist = cu.correlation_distribution( bc2_submat )

			#kdrew: temporarily add column to bicluster2, so we can calculate correlation, remove later
			bicluster2.add_column( i )
			bc2_submat = bicluster2.get_submatrix( self.__data_matrix )
			bc2_corrDist_with_i = cu.correlation_distribution( bc2_submat )
			bicluster2.remove_column( i )

			print "corrDist.mean(): %s" % (corrDist.mean(), )
			print "corrDist_without_i.mean(): %s" %(corrDist_without_i.mean(), )
			print "bc2_corrDist_with_i.mean(): %s" % (bc2_corrDist_with_i.mean(), )
			print "bc2_corrDist.mean(): %s" % (bc2_corrDist.mean(), )
			print ""

			corr_gain = corrDist.mean() - corrDist_without_i.mean()
			bc2_corr_gain = bc2_corrDist_with_i.mean() - bc2_corrDist.mean()

			print "corr_gain: %s" % (corr_gain,)
			print "bc2_corr_gain: %s" % (bc2_corr_gain,)
			
			return corr_gain, bc2_corr_gain



	#kdrew: when calling this on rows in bicluster, bicluster1 = bicluster, bicluster2 = random bicluster
	#kdrew: switch them when calling on rows from random bicluster
	def correlation_feature_row_by_index( self, bicluster1, bicluster2, index, i ):

			#kdrew: get portion of full data matrix that corresponds to the bicluster
			submat = bicluster1.get_submatrix( self.__data_matrix )

			#kdrew: calculate correlation of single row vs all other rows
			corrDistIndex = cu.correlation_distribution( submat, index )

			#kdrew: temporarily add row to random bicluster, so we can calculate correlation, remove later
			bicluster2.add_row(i)

			#kdrew: get portion of full matrix that corresponds to random bicluster
			bicluster2_submat = bicluster2.get_submatrix( self.__data_matrix )

			#kdrew: find index of row number in random bicluster
			#kdrew: ****this makes me a little nervous because we are getting the index from a set which is not ordered****
			bc2_index = list(bicluster2.rows()).index(i)

			#kdrew: calculate correlation of single row vs all other rows in random bicluster
			bc2_corrDistIndex =  cu.correlation_distribution( bicluster2_submat, bc2_index )

			bicluster2.remove_row(i)

			return corrDistIndex.mean(), bc2_corrDistIndex.mean()


