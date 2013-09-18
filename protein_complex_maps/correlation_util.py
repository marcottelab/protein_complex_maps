from scipy.stats.stats import pearsonr

def correlation_distribution(matrix, array_in=None):
	score_distribution = []
	pval_distribution = []

	for i, row in enumerate(matrix):
		#kdrew: if an array was passed in, compute correlation of all rows in matrix to array
		if array_in != None:
			score, pval = pearsonr(array_in, row)
			score_distribution.append(score)
			pval_distribution.append(pval)

		#kdrew: if no array is passed in, compute correlation of all rows in matrix to all other rows
		else:
			for j, row2 in enumerate(matrix[i+1:]):
				print i, row, j, row2
				score, pval = pearsonr(row, row2)
				score_distribution.append(score)
				pval_distribution.append(pval)

	return score_distribution, pval_distribution


