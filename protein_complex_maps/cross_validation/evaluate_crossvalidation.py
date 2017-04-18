
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import pickle
import pandas as pd
import itertools as it

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score

def main():

    parser = argparse.ArgumentParser(description="Evaluate Cross Validation, calc precision recall metrics")
    parser.add_argument("--results_wprob", action="store", dest="results_wprob", nargs='+', required=True, 
                                    help="Filenames of pair results with probability")
    parser.add_argument("--input_positives", action="store", dest="positives", nargs='+', required=True, 
                                    help="Filenames of positive pairs, files should be in order of results_wprob files")
    parser.add_argument("--input_negatives", action="store", dest="negatives", nargs='+', required=False, default=None,
                                    help="Filename of negative pairs, default = None (generated from processing positives), files should be in order of results_wprob files")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--threshold_recall", action="store", type=float, dest="threshold_recall", required=False, 
                                    help="Only calculate predictions at recall threshold")

    args = parser.parse_args()

    results_dict_list = []
    all_proteins_list = []
    ppi_list = []
    neg_ppi_list = []
    for i, filename in enumerate(args.results_wprob):
        results_dict = dict()
        results_file = open(filename,"rb")
        for line in results_file.readlines():
            if len(line.split()) == 3:
                id1 = line.split()[0]
                id2 = line.split()[1]
                ppi_str = str(sorted(list(frozenset([id1,id2]))))
                results_dict[ppi_str] = float(line.split()[2])
            else:
                print "WARNING: did not read line: %s" % (line,)

        #print "results size: %s" % (len(results_dict))
        results_dict_list.append(results_dict)


        positive_file = open(args.positives[i], "rb")
        all_proteins = set()
        ppis = set()
        neg_ppis = set()
        for line in positive_file.readlines():
            if len(line.split()) >= 2:
                id1 = line.split()[0]
                id2 = line.split()[1]
                all_proteins.add(id1)
                all_proteins.add(id2)

                ppi_str = str(sorted(list(frozenset([id1,id2]))))
                ppis.add(ppi_str)

        ppi_list.append(ppis)


        #kdrew: generate negative list by generating all pairs of proteins in positives but edges are not in positive list
        if args.negatives == None:
            for cpair in it.combinations(all_proteins,2):
                pair = str(sorted(list(frozenset(cpair))))
                if pair not in ppis:
                    #print "pair is neg: %s" % ' '.join(pair)
                    neg_ppis.add(pair)
        else:
            #kdrew: readin from file
            negative_file = open(args.negatives[i],"rb")
            for line in negative_file.readlines():
                if len(line.split()) >= 2:
                    id1 = line.split()[0]
                    id2 = line.split()[1]
                    ppi_str = str(sorted(list(frozenset([id1,id2]))))
                    neg_ppis.add(ppi_str)

        neg_ppi_list.append(neg_ppis)

    precision_list = []
    recall_list = []
    prauc_list = []
    for i, results_dict in enumerate(results_dict_list):
        true_array = []
        prob_array = []
        true_length = 0
        neg_length = 0
        for result_pair in results_dict.keys():
            if result_pair in ppi_list[i]:
                true_array.append(1)
                prob_array.append(results_dict[result_pair])
                true_length = true_length + 1

            elif result_pair in neg_ppi_list[i]:
                true_array.append(-1)
                prob_array.append(results_dict[result_pair])
                neg_length = neg_length + 1

        #print "Num Trues"
        #print str(true_length)
        #print "Num Negs"
        #print str(neg_length)
        #print len(true_array)
        #print len(prob_array)

        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        aps = average_precision_score(true_array, prob_array)
        prauc_list.append(aps)
        #print "precision: %s" % (list(precision),)
        #print "recall: %s" % (list(recall),)
        print "pr_auc: %s" % (aps,)

    print "mean pr auc: %s" % np.mean(prauc_list)


if __name__ == "__main__":
    main()


