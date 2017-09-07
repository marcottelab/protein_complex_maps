

import numpy as np
import pickle as p
import itertools as it

import pylab as pl
from scipy import interp

from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import StratifiedKFold
from sklearn.naive_bayes import GaussianNB

import protein_complex_maps.stoichiometry.relative_stoichiometry as rs


#kdrew: readin pdb benchmark for relative stoichiometries
stoich_file = open("/home/kdrew//scripts/protein_complex_maps/protein_complex_maps/stoichiometry/benchmark/ms_complete_pdbs.p")
stoich_data_set = p.load(stoich_file)

#kdrew: readin msds 
msds_file = open("/home/kdrew/data/protein_complex_maps/msds_pickles/HS_ms2_elutions_msds_peptide_normalized_ids_mapped.p")
#msds_file = open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/util/HS_ms2_elutions_msds_peptide_normalized_ids_mapped.p")
#msds_file = open("/home/kdrew/scripts/protein_complex_maps/protein_complex_maps/util/Hs_ms2_fromPeptide_msds_sc_mean_pdbs_ids_transfered.p")
msds = p.load(msds_file)


data_threshold = 50

all_stoichs = list()

data_list = list()
target_list = list()


#kdrew: setup data in n x 1 array, 1 feature of median of ratio list
for pdb in stoich_data_set:
	print pdb
	pdb_stoich = stoich_data_set[pdb]
	print pdb_stoich
	
	for pair in it.combinations(pdb_stoich.keys(),2):
		print pair
		pair_ratios = rs.calculate_ratio(msds, pair[0], pair[1])
		if len(pair_ratios) < data_threshold:
			continue
		median = np.median(pair_ratios)
		
		#kdrew: figure out real stoichiometry
		true_stoich = 1.0 * pdb_stoich[pair[0]] / pdb_stoich[pair[1]]

		try:
			stoich_index = all_stoichs.index(true_stoich)
		except ValueError:
			all_stoichs.append(true_stoich)
			stoich_index = all_stoichs.index(true_stoich)


		data_list.append(median)
		target_list.append(stoich_index)




for i in xrange(len(data_list)):
	print data_list[i], target_list[i]

data_array = np.array(data_list)
data_array = data_array.reshape(len(data_list),1)

target_array = np.array(target_list)

print all_stoichs


cv = StratifiedKFold(target_array, n_folds=6)
gnb = GaussianNB()
#gnb.fit(data_array, target_list)

mean_tpr = 0.0
mean_fpr = np.linspace(0,1,100)
all_tpr = []

for i, (train, test) in enumerate(cv):
	predictions = gnb.fit(data_array[train], target_array[train]).predict(data_array[test])

	print "test: %s" % target_array[test]
	print "predictions: %s" % predictions
	print "target_array: %s" % (target_array[test] == 0)
	print "predictions: %s" % (predictions == 0)
	fpr, tpr, thresholds = roc_curve(target_array[test], predictions == 0, pos_label=0)
	print "fpr: %s tpr: %s thresholds: %s" % (fpr,tpr, thresholds)
	mean_tpr += interp(mean_fpr, fpr, tpr)
	mean_tpr[0] = 0.0
	roc_auc = auc(fpr, tpr)
	pl.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

pl.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')

mean_tpr /= len(cv)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
pl.plot(mean_fpr, mean_tpr, 'k--',
label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

pl.xlim([-0.05, 1.05])
pl.ylim([-0.05, 1.05])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Receiver operating characteristic example')
pl.legend(loc="lower right")
pl.show()


#y_pred = gnb.predict(data_array)
#
#print "missed: %s out of %s" % ((target_list != y_pred).sum(), len(target_list))


