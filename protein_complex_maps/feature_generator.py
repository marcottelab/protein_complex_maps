
import numpy as np
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc
import pandas as pd
import statsmodels.api as sm

#kdrew: this class is for generating the features of bicluster rows and columns for creating the logistic regression model

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

	def create_logistic_regression(self, feature_matrix, training_cols):
		print feature_matrix
		print training_cols
		data = feature_matrix[training_cols]
		print "data: %s" % (data,)
		data['intercept'] = 1.0
		print "data: %s" % (data,)
		#kdrew: assumption label_header is 0, probably not(?)
		train_cols = data.columns[ data.columns != self.label_header ]
		print "train_cols: %s" % (train_cols,)
		logit = sm.Logit(data[self.label_header], data[train_cols])
		print logit
		result = logit.fit()
		print result
		return logit, result

	#kdrew: function to calculate the correlation feature for all columns in and out of bicluster, 
	#kdrew: creates a dataframe of features per bicluster column
	def correlation_feature_column(self, ):

		#kdrew: initialize new column and set all values to NAN
		self.__column_feature_matrix[self.corr_gain_bc_header] = float('nan')
		self.__column_feature_matrix[self.corr_gain_rand_header] = float('nan')

		#print self.__column_feature_matrix

		#kdrew: first do columns in bicluster and then do columns in random bicluster

		#kdrew: index is the index in the bicluster, i is the index (column number) in the complete data matrix
		for index, i in enumerate(self.__bicluster.columns()):

			#kdrew: think about removing a column from rand column bicluster to keep number of columns the same between rand and bc
			#print "index: %s, i: %s" % (index, i)
			corr_gain = self.correlation_feature_column_by_index( self.__bicluster, i )
			self.__rand_column_bicluster.add_column(i)
			rand_corr_gain = self.correlation_feature_column_by_index( self.__rand_column_bicluster, i )
			self.__rand_column_bicluster.remove_column(i)

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where column_header == i and column is corr_mean_ratio_header
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_bc_header ] = corr_gain
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_rand_header ] = rand_corr_gain

		for index, i in enumerate(self.__rand_column_bicluster.columns()):

			#kdrew: think about adding a column from rand column bicluster to keep number of columns the same between rand and bc
			rand_corr_gain = self.correlation_feature_column_by_index( self.__rand_column_bicluster, i )
			self.__bicluster.add_column(i)
			bc_corr_gain = self.correlation_feature_column_by_index( self.__bicluster, i )
			self.__bicluster.remove_column(i)

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all columns where column_header == i and column is corr_mean_ratio_header
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_bc_header ] = bc_corr_gain
			self.__column_feature_matrix.loc[ self.__column_feature_matrix[self.column_header]==i, self.corr_gain_rand_header ] = rand_corr_gain

		#print self.__column_feature_matrix


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

			corr_mean = self.correlation_feature_row_by_index( self.__bicluster, i )
			#kdrew: temporarily add row i, remove after calculation
			self.__rand_row_bicluster.add_row(i)
			rand_corr_mean = self.correlation_feature_row_by_index( self.__rand_row_bicluster, i )
			self.__rand_row_bicluster.remove_row(i)

			corr_mean_ratio = corr_mean/rand_corr_mean

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_bc_header ] = corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_rand_header ] = rand_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_ratio_header ] = corr_mean_ratio

		for index, i in enumerate(self.__rand_row_bicluster.rows()):

			rand_corr_mean = self.correlation_feature_row_by_index( self.__rand_row_bicluster, i )
			#kdrew: temporarily add row i, remove after calculation
			self.__bicluster.add_row(i)
			bc_corr_mean = self.correlation_feature_row_by_index( self.__bicluster, i )
			self.__bicluster.remove_row(i)

			rand_corr_mean_ratio = bc_corr_mean/rand_corr_mean

			#kdrew: this is kinda ugly but this is how you update values in a pandas dataframe
			#kdrew: find location (loc) by conditional: all rows where row_header == i and column is corr_mean_ratio_header
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_bc_header ] = bc_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_rand_header ] = rand_corr_mean
			self.__row_feature_matrix.loc[ self.__row_feature_matrix[self.row_header]==i, self.corr_mean_ratio_header ] = rand_corr_mean_ratio

		#print self.__row_feature_matrix

	#kdrew: when calling this on columns in bicluster, bicluster1 = bicluster, bicluster2 = random bicluster
	#kdrew: switch them when calling on columns from random bicluster
	def correlation_feature_column_by_index( self, bicluster1, i ):

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

			corr_gain = corrDist.mean() - corrDist_without_i.mean()

			return corr_gain


	#kdrew: when calling this on rows in bicluster, bicluster1 = bicluster, bicluster2 = random bicluster
	#kdrew: switch them when calling on rows from random bicluster
	def correlation_feature_row_by_index( self, bicluster1, i ):

			#kdrew: get portion of full data matrix that corresponds to the bicluster
			submat = bicluster1.get_submatrix( self.__data_matrix )

			index = list(bicluster1.rows()).index(i)

			#kdrew: calculate correlation of single row vs all other rows
			corrDistIndex = cu.correlation_distribution( submat, index )

			return corrDistIndex.mean()


