
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

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'grid': False})
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

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
    parser.add_argument("--labels", action="store", dest="labels", nargs='+', required=False, default=None,
                                    help="Labels for input results in order of --results_wprob list")
    parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, 
                                    help="Only tally predictions above probability threshold")
    parser.add_argument("--plot_thresholds", action="store_true", dest="plot_thresholds", required=False, default=False,
                                    help="Add probability threshold markers to plot")

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
        neg_list = []
        pos_list = []
        for result_pair in results_dict.keys():
            if args.threshold == None or args.threshold <= results_dict[result_pair]:
                if result_pair in ppis:
                    true_array.append(1)
                    prob_array.append(results_dict[result_pair])
                    pos_list.append((result_pair, results_dict[result_pair]))
                elif result_pair in neg_ppis:
                    true_array.append(-1)
                    prob_array.append(results_dict[result_pair])
                    neg_list.append((result_pair, results_dict[result_pair]))


        sorted_neg_list = sorted(neg_list, key=lambda k: k[1])
        sorted_pos_list = sorted(pos_list, key=lambda k: k[1])
        print "#ofNegs: %s, #ofPos: %s, FDR: %s" % (len(sorted_neg_list),len(sorted_pos_list),1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))
        for k in sorted_neg_list[-100:]:
            print "neg: %s:%s" % (k[0], k[1])
        for k in sorted_pos_list[-100:]:
            print "pos: %s:%s" % (k[0], k[1])
        print len(true_array)
        print len(prob_array)

        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        #average_precision = average_precision_score(true_array, prob_array)

        #print "precision: %s" % (list(precision),)
        #print "recall: %s" % (list(recall),)
        print len(precision)
        print len(recall)
        print len(thresholds)


        if args.plot_thresholds:
            #kdrew: set the index of the entry >= to threshold
            #kdrew: if there is no entry >= to threshold, set the index to the previous index
            try:
                indexOf1 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.1) 
            except StopIteration:
                pass
            try:
                indexOf2 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.2) 
            except StopIteration:
                indexOf2 = index01
            try:
                indexOf3 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.3) 
            except StopIteration:
                indexOf3 = indexOf2
            try:
                indexOf4 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.4) 
            except StopIteration:
                indexOf4 = indexOf3
            try:
                indexOf5 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.5) 
            except StopIteration:
                indexOf5 = indexOf4
            try:
                indexOf6 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.6) 
            except StopIteration:
                indexOf6 = indexOf5
            try:
                indexOf7 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.7) 
            except StopIteration:
                indexOf7 = indexOf6
            try:
                indexOf8 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.8) 
            except StopIteration:
                indexOf8 = indexOf7
            try:
                indexOf9 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.9) 
            except StopIteration:
                indexOf9 = indexOf8
            try:
                indexOf10 = next(x[0] for x in enumerate(thresholds) if x[1] >= 1.0) 
            except StopIteration:
                indexOf10 = indexOf9

            print thresholds[indexOf1]
            print indexOf1
            print thresholds[indexOf2]
            print indexOf2
            print thresholds[indexOf3]
            print indexOf3
            print thresholds[indexOf4]
            print indexOf4
            print thresholds[indexOf10]
            print indexOf10

            #print thresholds[indexOf1-1]
            #print thresholds[indexOf10-1]

            threshold_precisions = []
            threshold_recalls = []
            threshold_labels = []
            print "threshold 0.1"
            print "precision: %s" % precision[indexOf1]
            print "recall: %s" % recall[indexOf1]
            threshold_precisions.append(precision[indexOf1])
            threshold_recalls.append(recall[indexOf1])
            threshold_labels.append('0.1')
            
            print "threshold 0.2"
            print "precision: %s" % precision[indexOf2]
            print "recall: %s" % recall[indexOf2]
            threshold_precisions.append(precision[indexOf2])
            threshold_recalls.append(recall[indexOf2])
            threshold_labels.append('0.2')
            
            print "threshold 0.3"
            print "precision: %s" % precision[indexOf3]
            print "recall: %s" % recall[indexOf3]
            threshold_precisions.append(precision[indexOf3])
            threshold_recalls.append(recall[indexOf3])
            threshold_labels.append('0.3')
            
            print "threshold 0.4"
            print "precision: %s" % precision[indexOf4]
            print "recall: %s" % recall[indexOf4]
            threshold_precisions.append(precision[indexOf4])
            threshold_recalls.append(recall[indexOf4])
            threshold_labels.append('0.4')
            
            print "threshold 0.5"
            print "precision: %s" % precision[indexOf5]
            print "recall: %s" % recall[indexOf5]
            threshold_precisions.append(precision[indexOf5])
            threshold_recalls.append(recall[indexOf5])
            threshold_labels.append('0.5')
            
            print "threshold 0.6"
            print "precision: %s" % precision[indexOf6]
            print "recall: %s" % recall[indexOf6]
            threshold_precisions.append(precision[indexOf6])
            threshold_recalls.append(recall[indexOf6])
            threshold_labels.append('0.6')
            
            print "threshold 0.7"
            print "precision: %s" % precision[indexOf7]
            print "recall: %s" % recall[indexOf7]
            threshold_precisions.append(precision[indexOf7])
            threshold_recalls.append(recall[indexOf7])
            threshold_labels.append('0.7')

            print "threshold 0.8"
            print "precision: %s" % precision[indexOf8]
            print "recall: %s" % recall[indexOf8]
            threshold_precisions.append(precision[indexOf8])
            threshold_recalls.append(recall[indexOf8])
            threshold_labels.append('0.8')
            
            print "threshold 0.9"
            print "precision: %s" % precision[indexOf9]
            print "recall: %s" % recall[indexOf9]
            threshold_precisions.append(precision[indexOf9])
            threshold_recalls.append(recall[indexOf9])
            threshold_labels.append('0.9')

            print "threshold 1.0"
            print "precision: %s" % precision[indexOf10]
            print "recall: %s" % recall[indexOf10]
            threshold_precisions.append(precision[indexOf10])
            threshold_recalls.append(recall[indexOf10])
            threshold_labels.append('1.0')


        label = args.results_wprob[i]
        if args.labels != None:
            label = args.labels[i]
        line, = plt.plot(recall, precision, label=label)
        if args.plot_thresholds:
            plt.scatter(threshold_recalls, threshold_precisions, color=line.get_color())
            for i, label in enumerate(threshold_labels):
                plt.annotate(label, xy=(threshold_recalls[i],threshold_precisions[i]))


    #plt.clf()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    #plt.title('Precision-Recall example: AUC={0:0.2f}'.format(average_precision))
    plt.title('Precision-Recall')
    plt.legend(loc="upper right",fontsize=8)

    plt.savefig(args.output_file)


if __name__ == "__main__":
    main()


