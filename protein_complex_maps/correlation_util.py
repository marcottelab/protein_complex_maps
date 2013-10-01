from scipy.stats.stats import pearsonr
import numpy as np

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


