from scipy.stats.stats import pearsonr
import math as m
import numpy as np

CLOSE_TO_ONE = 0.99999999
EPSILON = 0.000001

def correlation_distribution(matrix, index=None, nans=0):

	corrcoefMat = np.corrcoef(matrix)
	#print corrcoefMat
	if index == None:
		#kdrew: gets indices of upper triangle (w/o diagonal) and returns values in a list
		#score_distribution = corrcoefMat[np.triu_indices(len(corrcoefMat),1)].tolist()
		score_distribution = corrcoefMat[np.triu_indices(len(corrcoefMat),1)]

	else:
		#kdrew: return the correlation scores of the row identified by index to all other rows, do not include to itself
		score_distribution = np.delete(corrcoefMat[index],index)

	#print score_distribution

	if nans != None:
		nan_list = np.isnan(score_distribution)
		#print nan_list
		try:
			score_distribution[nan_list] = nans
		except TypeError:
			pass

	return score_distribution


#kdrew: calculates t-values for array of correlation coefficents
def tvalue_correlation(ar, n):

	#kdrew: t-value formula
	def tvalue(r,n):
		try:
			#kdrew: get a divide by zero if r is 1.0 (perfect correlation), set r to be close to one to avoid that
			if m.fabs(r-1.0) < EPSILON: 
				r = CLOSE_TO_ONE
				print "tvalue_correlation: r == 1.0, changing r value to be CLOSE_TO_ONE" 
			if n <= 2:
				print "tvalue_correlation: degrees of freedom <= 0, returning tvalue = 0.0"
				return 0.0

			t = r/m.sqrt((1-r**2)/(n-2))

		except ZeroDivisionError:
			print "ZeroDivisionError: r: %s n %s" % (r, n)
			t = 0.0
		except ValueError:
			print "ValueError: r: %s n %s" % (r, n)
			t = 0.0


		print "r: %s n: %s t: %s" % (r, n, t)
		return t
		
	tvalue_vectorized = np.vectorize(tvalue)
	return tvalue_vectorized(ar, n)

