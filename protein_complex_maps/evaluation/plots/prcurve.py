
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
from sklearn.metrics import auc


from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
#rcParams.update({'grid': False})
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
    parser.add_argument("--threshold_markers", action="store", dest="threshold_markers", nargs='+', required=False, default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],
                                    help="Sets which thresholds to plot along prcurve, default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]")
    parser.add_argument("--precision_markers", action="store", dest="precision_markers", nargs='+', required=False, default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],
                                    help="Sets which precisions to report recall, default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]")
    parser.add_argument("--complete_benchmark", action="store_true", dest="complete_benchmark", required=False, default=False,
                                    help="Use the complete benchmark and set the probablility to 0.0, default=False")
    parser.add_argument("--add_tiny_noise", action="store_true", dest="add_tiny_noise", required=False, default=False,
                                    help="Add tiny bit of noise to scores of each prediction to separate predictions with same score, default=False")

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

        label = args.results_wprob[i]
        if args.labels != None:
            label = args.labels[i]
        print("label: %s" % label)

        for result_pair in results_dict.keys():
            #kdrew: add a tiny bit of noise to scores so that predictions with all of the same scores are separated, should fix straight line at end of pr plots
            if args.add_tiny_noise:
                tiny_noise = np.random.random()/1000000
            else:
                tiny_noise = 0.0
            if args.threshold == None or args.threshold <= results_dict[result_pair]:
                if result_pair in ppis:
                    true_array.append(1)
                    prob_array.append(results_dict[result_pair]+tiny_noise)
                    pos_list.append((result_pair, results_dict[result_pair]))
                elif result_pair in neg_ppis:
                    true_array.append(-1)
                    prob_array.append(results_dict[result_pair]+tiny_noise)
                    neg_list.append((result_pair, results_dict[result_pair]))

        if args.complete_benchmark:
            for ppi in set(ppis) - set(results_dict.keys()):
                print "complete_benchmark pos: %s" % ppi
                true_array.append(1)
                if args.add_tiny_noise:
                    tiny_noise = np.random.random()/1000000
                    prob_array.append(tiny_noise)
                else:
                    prob_array.append(0.0)
            for neg_ppi in set(neg_ppis) - set(results_dict.keys()):
                true_array.append(-1)
                if args.add_tiny_noise:
                    tiny_noise = np.random.random()/1000000
                    prob_array.append(tiny_noise)
                else:
                    prob_array.append(0.0)


        sorted_neg_list = sorted(neg_list, key=lambda k: k[1])
        sorted_pos_list = sorted(pos_list, key=lambda k: k[1])
        print "#ofNegs: %s" % (len(sorted_neg_list))
        print "#ofPos: %s" % (len(sorted_pos_list))
        print "FDR: %s" % (1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))

        print "#ofNegs: %s, #ofPos: %s, FDR: %s" % (len(sorted_neg_list),len(sorted_pos_list),1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))
        for k in sorted_neg_list[-100:]:
            print "neg: %s:%s" % (k[0], k[1])
        for k in sorted_pos_list[-100:]:
            print "pos: %s:%s" % (k[0], k[1])
        print len(true_array)
        print len(prob_array)

        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        average_precision = average_precision_score(true_array, prob_array)
        area_under_curve = auc(recall, precision)

        #print "precision: %s" % (list(precision),)
        #print "recall: %s" % (list(recall),)
        print len(precision)
        print len(recall)
        print len(thresholds)
        print("average_precision: %s" % average_precision)
        print("area_under_curve: %s" % area_under_curve)

        #kdrew: from https://www.geeksforgeeks.org/python-find-closest-number-to-k-in-given-list/ 

        for x in args.precision_markers:
            closest_prec = precision[min(range(len(precision)), key = lambda i: abs(precision[i]-float(x)))]
            index_of_prec = list(precision).index(closest_prec)
            try:
                print "precision: %s, closest_precision: %s, recall: %s, threshold: %s" % (x, closest_prec, recall[index_of_prec], thresholds[index_of_prec])
            except IndexError, e:
                print e
                continue



        if args.plot_thresholds:
            #kdrew: set the index of the entry >= to threshold
            threshold_indices = dict()
            #kdrew: ensure passed in arguments are floats
            for threshold_marker in [float(x) for x in args.threshold_markers]:
                try:
                    threshold_indices[threshold_marker] = next(i for i,x in enumerate(thresholds) if x >= threshold_marker) 
                except StopIteration:
                    pass

            threshold_precisions = []
            threshold_recalls = []
            threshold_labels = []
            for threshold_marker in threshold_indices:
                threshold_precisions.append(precision[threshold_indices[threshold_marker]])
                threshold_recalls.append(recall[threshold_indices[threshold_marker]])
                threshold_labels.append(threshold_marker)

                print "threshold %s" % threshold_marker
                print "precision: %s" % precision[threshold_indices[threshold_marker]]
                print "recall: %s" % recall[threshold_indices[threshold_marker]]


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


