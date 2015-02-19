
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

parser = argparse.ArgumentParser(description="")

parser.add_argument("--results_file", action="store", dest="results_filename", required=True, help="")

args = parser.parse_args()

f = open(args.results_filename, "rb")

bc_zscores = []
interaction_zscores = []
for line in f.readlines():
	bc_zscore = float(line.split()[2])
	int_zscore = float(line.split("interaction_zscore:")[1])

	if np.isnan(int_zscore) or np.isnan(bc_zscore):
		continue
	bc_zscores.append(bc_zscore)
	interaction_zscores.append(int_zscore)

print "pearsonr: %s" % (pearsonr(bc_zscores, interaction_zscores),)
print "spearmanr: %s" % (spearmanr(bc_zscores, interaction_zscores),)
plt.scatter(bc_zscores, interaction_zscores)
plt.show()

