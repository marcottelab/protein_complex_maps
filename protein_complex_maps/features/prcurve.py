
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

    parser = argparse.ArgumentParser(description="Generate Precision Recall Curve")
    parser.add_argument("--results_wprob", action="store", dest="results_wprob", nargs='+', required=True, 
                                    help="Filename of pair results with probability")
    parser.add_argument("--input_positives", action="store", dest="positives", required=True, 
                                    help="Filename of positive pairs")
    parser.add_argument("--input_negatives", action="store", dest="negatives", required=False, default=None,
                                    help="Filename of negative pairs, default = None (generated from processing positives)")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")

    args = parser.parse_args()

    results_dict_list = []
    for filename in args.results_wprob:
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

        print "results size: %s" % (len(results_dict))
        results_dict_list.append(results_dict)

    positive_file = open(args.positives,"rb")

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


    #kdrew: generate negative list by generating all pairs of proteins in positives but edges are not in positive list
    if args.negatives == None:
        for cpair in it.combinations(all_proteins,2):
            pair = str(sorted(list(frozenset(cpair))))
            if pair not in ppis:
                #print "pair is neg: %s" % ' '.join(pair)
                neg_ppis.add(pair)
    else:
        #kdrew: readin from file
        negative_file = open(args.negatives,"rb")
        for line in negative_file.readlines():
            if len(line.split()) >= 2:
                id1 = line.split()[0]
                id2 = line.split()[1]
                ppi_str = str(sorted(list(frozenset([id1,id2]))))
                neg_ppis.add(ppi_str)


    for i, results_dict in enumerate(results_dict_list):
        true_array = []
        prob_array = []
        for result_pair in results_dict.keys():
            if result_pair in ppis:
                true_array.append(1)
                prob_array.append(results_dict[result_pair])
            elif result_pair in neg_ppis:
                true_array.append(-1)
                prob_array.append(results_dict[result_pair])


        print len(true_array)
        print len(prob_array)

        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        #average_precision = average_precision_score(true_array, prob_array)

        plt.plot(recall, precision, label=args.results_wprob[i])


    #plt.clf()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    #plt.title('Precision-Recall example: AUC={0:0.2f}'.format(average_precision))
    plt.title('Precision-Recall')
    plt.legend(loc="lower left",fontsize=8)

    plt.savefig(args.output_file)


if __name__ == "__main__":
    main()


