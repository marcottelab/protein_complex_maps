
import sys
import numpy as np
import random as r
import protein_complex_maps.correlation_util as cu
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.feature_generator as fg
from scipy.stats import ks_2samp
import statsmodels.tools.sm_exceptions
import protein_complex_maps.plots.plot_bicluster as pb
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.score_util as su
import protein_complex_maps.read_data as rd
import protein_complex_maps.bicluster_generator as bg


#kdrew: there are several ways to massage these data
# remove columns with all zeros
# add pinch of noise (1/M) 
# normalize so that whole column sums to one (this essentially makes all columns equal to each other)
# scale values * 1/(min value) so all values are >= 1 (helps with multiplication scores that converge to single row, single column, scale after normalize that makes all values <= 1)
# apply poisson distribution (this puts more variance to larger values and therefore makes larger values noiser than smaller values)
# apply gaussian distribution (this puts equal variance to all values, need to choose sigma (default 1)
# applying distributions can be coupled with multiple rounds of sampling and scoring, the average of scores is returned



r.seed(123456)
np.random.seed(1234)
sample_filename1 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"
sample_filename2 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome_exosc_randos.txt"

sample_iterations = 100

##kdrew: this may be unnecessary because the default sigma is 1.0
#def normal_sigma1(mu):
#   return np.random.normal(mu,1.0)
	

sample_file1 = open(sample_filename1, 'rb')
sample_file2 = open(sample_filename2, 'rb')

data_matrix1, name_list1 = rd.read_datafile(sample_file1)
data_matrix2, name_list2 = rd.read_datafile(sample_file2)

##kdrew: eat header
#line = sample_file.readline()
#
#data = []
#
#for line in sample_file.readlines():
#	#print line
#	line_data = line.split()
#	line_array = map(float,line_data[2:])
#	print line_array
#	data.append(line_array)
#
##print data
#
#data_matrix = np.asmatrix(data)
#
##print data_matrix


#kdrew: remove columns with all zeros
#clean_data_matrix = data_matrix.compress(~np.array(np.all(data_matrix[:]==0,axis=0))[0],axis=1)
clean_data_matrix_pre_normalized1 = nu.remove_zero(data_matrix1)
clean_data_matrix_pre_normalized2 = nu.remove_zero(data_matrix2)
#print clean_data_matrix_pre_normalized1
##kdrew: remove columns with all zeros
#self.__clean_data_matrix = self.__data_matrix.compress(~np.array(np.all(self.__data_matrix[:]==0,axis=0))[0],axis=1)

#print "percentile1: %s" % (np.percentile(clean_data_matrix_pre_normalized1,75),)
#print "percentile2: %s" % (np.percentile(clean_data_matrix_pre_normalized2,75),)

#clean_data_matrix_binary1 = nu.binary(clean_data_matrix_pre_normalized1, np.percentile(clean_data_matrix_pre_normalized1,75) )
#clean_data_matrix_binary2 = nu.binary(clean_data_matrix_pre_normalized2, np.percentile(clean_data_matrix_pre_normalized2,75) )

clean_data_matrix_binary1 = nu.binary(clean_data_matrix_pre_normalized1, 1 )
clean_data_matrix_binary2 = nu.binary(clean_data_matrix_pre_normalized2, 1 )

#clean_data_matrix_noised1 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized1)
#clean_data_matrix_noised2 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized2)
clean_data_matrix_noised1 = nu.add_noise_over_columns(clean_data_matrix_binary1)
clean_data_matrix_noised2 = nu.add_noise_over_columns(clean_data_matrix_binary2)

#kdrew: normalize where whole column adds to 1.0
#clean_data_matrix = nu.normalize_over_rows(clean_data_matrix_noised)
clean_data_matrix, name_list = rd.concat_data_matrix( clean_data_matrix_noised1, name_list1, clean_data_matrix_noised2, name_list2)
#clean_data_matrix, name_list = rd.concat_data_matrix( clean_data_matrix_binary1, name_list1, clean_data_matrix_binary2, name_list2)
print clean_data_matrix

#clean_data_matrix = nu.min_to_one_scale(clean_data_matrix)
#print clean_data_matrix

#bicluster1 = bc.Bicluster(rows = [3,10,13,14,22,26], cols = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70], random_module=r)
#bicluster1 = bc.Bicluster(rows = [3,9], cols = [20,22], random_module=r)

#bc1_submat = bicluster1.get_submatrix(clean_data_matrix)
#print bc1_submat

#bc1_corrDist = cu.correlation_distribution(bc1_submat)
#print bc1_corrDist
#
#bc1_r0_corrDist = cu.correlation_distribution(bc1_submat, 0)
#print bc1_r0_corrDist
#
#rand_bicluster1 = bc.Bicluster(rows=bicluster1.get_random_outside_rows(clean_data_matrix), cols=bicluster1.columns(), random_module=r)
#
#kdrew: index is the index in the bicluster, i is the index (row number) in the complete data matrix
#for index, i in enumerate(bicluster1.rows()):
#	print "index:%s, i:%s" % (index, i)
#	print bicluster1.get_row(clean_data_matrix, row=i)
#	submat = bicluster1.get_submatrix( clean_data_matrix )
#	print "bicluster1 submatrix:"
#	print submat
#
#	print ""
#	submatWO = bicluster1.get_submatrix( clean_data_matrix, without_rows=[i] )
#	corrDistWO = cu.correlation_distribution( submatWO )
#	print "correlation_distribution without %s:" % (i,)
#	print corrDistWO
#
#	corrDistIndex = cu.correlation_distribution( submat, index )
#	print "correlation_distribution of bicluster submatrix vs index %s:" % (index,)
#	print corrDistIndex
#
#	#ks_biclust = ks_2samp(corrDistWO, corrDistIndex)
#	#print ks_biclust
#
#	rand_bicluster1.add_row(i)
#
#	rand_submat = rand_bicluster1.get_submatrix( clean_data_matrix )
#	print "random bicluster submatrix with row %s:" % (i,)
#	print rand_submat
#	rand_submatWO = rand_bicluster1.get_submatrix( clean_data_matrix, without_rows=[i] )
#
#	randbc_corrDistWO = cu.correlation_distribution( rand_submatWO )
#	print randbc_corrDistWO
#
#	rand_bc_index = list(rand_bicluster1.rows()).index(i)
#	randbc_corrDistIndex =  cu.correlation_distribution( rand_submat, rand_bc_index )
#	print randbc_corrDistIndex
#
#	#kdrew: not really using ks tests anymore, 10/09/13
#	#ks_rand_biclust = ks_2samp(randbc_corrDistWO, randbc_corrDistIndex)
#	#print ks_rand_biclust
#
#	#kdrew: compute ratio of ks tests, this gives a value of how correlated the bicluster is compared to the background
#	#kdrew: values of ~1.0 suggest no correlation above background, small values suggest bicluster is more correlated than background
#	#print ks_biclust[0]/ks_rand_biclust[0]
#
#	print corrDistIndex.mean()
#	print randbc_corrDistIndex.mean()
#	print corrDistIndex.mean()/randbc_corrDistIndex.mean()
#
#	print ""
#
#	rand_bicluster1.remove_row(i)


#featuregen = fg.FeatureGenerator(bicluster1, clean_data_matrix)
#featuregen.correlation_feature_row()
#featuregen.correlation_feature_column()
#
#test_columns = ['label','corr_mean_bc','corr_mean_rand']
#try:
#	row_logreg, row_logit_result = featuregen.create_logistic_regression_row(test_columns)
#	print row_logit_result.summary()
#except statsmodels.tools.sm_exceptions.PerfectSeparationError:
#	print "PerfectSeparationError"
#	sys.exit(0)
#
#try:
#	column_logreg, column_logit_result = featuregen.create_logistic_regression_column()
#	print column_logit_result.summary()
#except statsmodels.tools.sm_exceptions.PerfectSeparationError:
#	print "PerfectSeparationError"
#	sys.exit(0)
# 

##kdrew: annealing temperature
#T=4
##kdrew: randomly pick row or column (inside or outside of bicluster)
#print "rows: %s" % (len(clean_data_matrix),)
#numRows, numCols = clean_data_matrix.shape
#print numRows, numCols
#
##kdrew: probably should test for some convergence
#for i in xrange(1,1000):
#	print "iteration: %s" % (i,)
#
#	random_row = np.random.random_integers(0, numRows-1)
#	random_column = np.random.random_integers(0, numCols-1)
#
#	print "bicluster matrix: %s" % (bicluster1.get_submatrix(clean_data_matrix))
#
##	featuregen = fg.FeatureGenerator(bicluster1, clean_data_matrix)
##	featuregen.correlation_feature_row()
##	featuregen.correlation_feature_column()
##
##	#test_columns = ['label','corr_mean_bc','corr_mean_rand']
##	#test_columns = ['label','corr_mean_bc','tvalue_corr_mean_bc']
##	test_columns = ['label','tvalue_corr_mean_bc']
##
##	try:
##		row_logreg, row_logit_result = featuregen.create_logistic_regression_row(test_columns)
##		try:
##			print row_logit_result.summary()
##		except ValueError:
##			pass
##	except statsmodels.tools.sm_exceptions.PerfectSeparationError:
##		print "PerfectSeparationError"
##
##	#test_columns = ['label','corr_gain_bc','tvalue_corr_gain_bc']
##	test_columns = ['label','tvalue_corr_gain_bc']
##	try:
##		column_logreg, column_logit_result = featuregen.create_logistic_regression_column(test_columns)
##		try:
##			print column_logit_result.summary()
##		except ValueError:
##			pass
##	except statsmodels.tools.sm_exceptions.PerfectSeparationError:
##		print "PerfectSeparationError"
##
#
#	if random_row in bicluster1.rows():
#		print "row %s in bicluster" % (random_row,)
#		##kdrew: if row is inside bicluster, test if it decreases total mean correlation
#		#corr_with_row = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_with_row = corr_with_row.mean()
#		##kdrew: calculate tvalues for correlation coefficients
#		#tvalue_corr_with_row = cu.tvalue_correlation(corr_with_row, len(bicluster1.columns()))
#		##mean_tvalue_corr_with_row = tvalue_corr_with_row.mean()
#
#		#normal_corr_with_row, normal_tvalue_corr_with_row = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_with_row = normal_tvalue_corr_with_row.mean()
#
#		multidot_score_with_row = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_with_row = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.remove_row(random_row)
#		print "bicluster matrix without row: %s" % (bicluster1.get_submatrix(clean_data_matrix))
#
#
#		#corr_without_row = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_without_row = corr_without_row.mean()
#		#tvalue_corr_without_row = cu.tvalue_correlation(corr_without_row, len(bicluster1.columns()))
#		##mean_tvalue_corr_without_row = tvalue_corr_without_row.mean()
#
#		#normal_corr_without_row, normal_tvalue_corr_without_row = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_without_row = normal_tvalue_corr_without_row.mean()
#
#		multidot_score_without_row = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_without_row = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#
#		bicluster1.add_row(random_row)
#
#		#kdrew: I was first using mean correlation as a metric to add or remove a row, now using tvalues
#		#if(mean_corr_with_row >= mean_corr_without_row):
#		#if(mean_tvalue_corr_with_row >= mean_tvalue_corr_without_row):
#		print "score with: %s without: %s" % (multidot_score_with_row, multidot_score_without_row)
#		if( multidot_score_with_row >= multidot_score_without_row ):
#			print "row does not decrease total correlation, test to keep"
#			#kdrew: if no, use logistic regression to predict membership value: x
#			#kdrew: p(drop|x) = math.e**(-(1-x)/T)
#
#		else:
#			#kdrew: if yes, remove row #			print "row decreases total correlation, automatically remove"
#			bicluster1.remove_row(random_row)
#
#	else:
#		#kdrew: if row is outside bicluster, test if it increases total mean correlation, 
#		print "row %s outside of bicluster" % (random_row,)
#		#corr_without_row = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_without_row = corr_without_row.mean()
#		##kdrew: calculate tvalues for correlation coefficients
#		#tvalue_corr_without_row = cu.tvalue_correlation(corr_without_row, len(bicluster1.columns()))
#		##mean_tvalue_corr_without_row = tvalue_corr_without_row.mean()
#
#		#normal_corr_without_row, normal_tvalue_corr_without_row = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_without_row = normal_tvalue_corr_without_row.mean()
#
#		multidot_score_without_row = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_without_row = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.add_row(random_row)
#		print "bicluster matrix with row: %s" % (bicluster1.get_submatrix(clean_data_matrix))
#
#		#corr_with_row = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_with_row = corr_with_row.mean()
#		##kdrew: calculate tvalues for correlation coefficients
#		#tvalue_corr_with_row = cu.tvalue_correlation(corr_with_row, len(bicluster1.columns()))
#		##mean_tvalue_corr_with_row = tvalue_corr_with_row.mean()
#
#		#normal_corr_with_row, normal_tvalue_corr_with_row = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_with_row = normal_tvalue_corr_with_row.mean()
#
#		multidot_score_with_row = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_with_row = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.remove_row(random_row)
#
#
#		#print mean_corr_with_row, mean_corr_without_row
#
#		#kdrew: I was first using mean correlation as a metric to add or remove a row, now using tvalues
#		#if(mean_corr_with_row > mean_corr_without_row):
#		#if(mean_tvalue_corr_with_row > mean_tvalue_corr_without_row):
#		print "score with: %s without: %s" % (multidot_score_with_row, multidot_score_without_row)
#		if(multidot_score_with_row > multidot_score_without_row):
#			#kdrew: if yes, add row
#			print "row increases total correlation, automatically add"
#			bicluster1.add_row(random_row)
#		else:
#			print "row does not increase total correlation, test to add"
#			#kdrew: if no, use logistic regression to predict membership value: x
#			#kdrew: p(add|x) = math.e**(-x/T)
#
#	#kdrew: check to make sure random column is not all zeros across the rows of the bicluster, if so remove or do not add
#	if np.all(bicluster1.get_column(clean_data_matrix, random_column) == 0):
#		print "removing: random column is all zeros: %s" % (bicluster1.get_column(clean_data_matrix,random_column),)
#		bicluster1.remove_column(random_column)
#		continue
#
#	
#	if random_column in bicluster1.columns():
#		print "column %s in bicluster" % (random_column,)
#		##kdrew: if column is inside bicluster, test if it decreases total mean correlation
#		#corr_with_column = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_with_column = corr_with_column.mean()
#		#tvalue_corr_with_column = cu.tvalue_correlation(corr_with_column, len(bicluster1.columns()))
#		##mean_tvalue_corr_with_column = tvalue_corr_with_column.mean()
#
#		#normal_corr_with_column, normal_tvalue_corr_with_column = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_with_column = normal_tvalue_corr_with_column.mean()
#
#		multidot_score_with_column = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_with_column = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.remove_column(random_column)
#		print "bicluster matrix without column: %s" % (bicluster1.get_submatrix(clean_data_matrix))
#
#		#corr_without_column = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_without_column = corr_without_column.mean()
#		#tvalue_corr_without_column = cu.tvalue_correlation(corr_without_column, len(bicluster1.columns()))
#		##mean_tvalue_corr_without_column = tvalue_corr_without_column.mean()
#
#		#normal_corr_without_column, normal_tvalue_corr_without_column = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_without_column = normal_tvalue_corr_without_column.mean()
#
#		multidot_score_without_column = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_without_column = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.add_column(random_column)
#
#		#if(mean_corr_with_column >= mean_corr_without_column ):
#		#if(mean_tvalue_corr_with_column >= mean_tvalue_corr_without_column ):
#		print "score with: %s without: %s" % (multidot_score_with_column, multidot_score_without_column)
#		if(multidot_score_with_column >= multidot_score_without_column):
#			print "column does not decrease total correlation, test to keep"
#			#kdrew: if no, use logistic regression to predict membership value: x
#			#kdrew: p(drop|x) = math.e**(-(1-x)/T)
#
#		else:
#			#kdrew: if yes, remove column
#			print "column decreases total correlation, automatically remove"
#			bicluster1.remove_column(random_column)
#
#	else:
#		#kdrew: if column is outside bicluster, test if it increases total mean correlation, 
#		print "column %s outside of bicluster" % (random_column,)
#		#corr_without_column = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_without_column = corr_without_column.mean()
#		#tvalue_corr_without_column = cu.tvalue_correlation(corr_without_column, len(bicluster1.columns()))
#		##mean_tvalue_corr_without_column = tvalue_corr_without_column.mean()
#
#		#normal_corr_without_column, normal_tvalue_corr_without_column = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_without_column = normal_tvalue_corr_without_column.mean()
#
#		multidot_score_without_column = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_without_column = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.add_column(random_column)
#		print "bicluster matrix with column: %s" % (bicluster1.get_submatrix(clean_data_matrix))
#
#		#corr_with_column = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#		#mean_corr_with_column = corr_with_column.mean()
#		#tvalue_corr_with_column = cu.tvalue_correlation(corr_with_column, len(bicluster1.columns()))
#		##mean_tvalue_corr_with_column = tvalue_corr_with_column.mean()
#
#		#normal_corr_with_column, normal_tvalue_corr_with_column = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#		#mean_tvalue_corr_with_column = normal_tvalue_corr_with_column.mean()
#
#		multidot_score_with_column = su.multiple_dot(bicluster1.get_submatrix(clean_data_matrix))
#		#multidot_score_with_column = su.sum_cells(bicluster1.get_submatrix(clean_data_matrix))
#		bicluster1.remove_column(random_column)
#
#		#print mean_corr_with_column, mean_corr_without_column
#		#if(mean_corr_with_column > mean_corr_without_column):
#		#if(mean_tvalue_corr_with_column > mean_tvalue_corr_without_column):
#		print "score with: %s without: %s" % (multidot_score_with_column, multidot_score_without_column)
#		if(multidot_score_with_column > multidot_score_without_column):
#			#kdrew: if yes, add column
#			print "column increases total correlation, automatically add"
#			bicluster1.add_column(random_column)
#		else:
#			print "column does not increase total correlation, test to add"
#			#kdrew: if no, use logistic regression to predict membership value: x
#			#kdrew: p(add|x) = math.e**(-x/T)
#


#kdrew: scale cells of bicluster in full matrix so we don't repeatidly generate the same biclusters
scale_factor = 0.75
working_data_matrix = clean_data_matrix

bcgen = bg.BiclusterGenerator(su.multiple_dot_neg, iterations=2500, random_module=np.random)
for i in xrange(1,100):
	bicluster1 = bcgen.generator(working_data_matrix)

	print "bicluster%s" % (i,)
	print bicluster1.rows()
	print bicluster1.columns()
	#print bicluster1.get_submatrix(working_data_matrix)

	#kdrew: not sure if evaluating on the clean_data_matrix is okay because the bicluster was optimized on the working_data_matrix
	eval_dict = bcgen.evaluate( clean_data_matrix, len(bcgen.biclusters)-1 )
	for t in eval_dict.keys():
		print "random %s, mean: %s, std: %s, zscore: %s" % ( t, eval_dict[t]['mean'], eval_dict[t]['std'], eval_dict[t]['zscore'] )

	pb.plot_bicluster(clean_data_matrix, bicluster1, savefilename="/home/kdrew/public_html/test/bicluster%s_plot.pdf" % (i))

	#kdrew: this value is a little arbitrary but only scale "converged" biclusters, there probably is a better way to set this or evaluate this
	if eval_dict['all']['zscore'] > 100:
		print "bicluster converged, scaling"
		working_data_matrix = bicluster1.scale(working_data_matrix, scale_factor)

#bc_normal_corr, bc_normal_tvalue_corr = cu.sample_correlation_distribution(bicluster1.get_submatrix(clean_data_matrix), noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#print "bicluster normal_corr.mean: %s normal_tvalue_corr.mean: %s" % (bc_normal_corr.mean(), bc_normal_tvalue_corr.mean())
#
#whole_normal_corr, whole_normal_tvalue_corr = cu.sample_correlation_distribution(clean_data_matrix, noise_constant=1.0/clean_data_matrix.shape[1], sample_module = normal_sigma1, iterations=sample_iterations)
#print "whole matrix normal_corr.mean: %s normal_tvalue_corr.mean: %s" % (whole_normal_corr.mean(), whole_normal_tvalue_corr.mean())
#
#bc_corr = cu.correlation_distribution(bicluster1.get_submatrix(clean_data_matrix))
#bc_tvalue = cu.tvalue_correlation(bc_corr, len(bicluster1.columns()))
#print "bicluster matrix corr.mean: %s tvalue_corr.mean: %s" % (bc_corr.mean(), bc_tvalue.mean())
#
#whole_corr = cu.correlation_distribution(clean_data_matrix)
#whole_tvalue = cu.tvalue_correlation(whole_corr, clean_data_matrix.shape[1])
#print "whole matrix corr.mean: %s tvalue_corr.mean: %s" % (whole_corr.mean(), whole_tvalue.mean())





