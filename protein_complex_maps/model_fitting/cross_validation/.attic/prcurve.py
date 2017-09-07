
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


    #tmp_outfilepos = open('atobsc3_arathtraesoryshbraoselml_pos.txt', "w")
    #tmp_outfileneg = open('atobsc3_arathtraesoryshbraoselml_neg.txt', "w")

    for i, results_dict in enumerate(results_dict_list):
        true_array = []
        prob_array = []
        true_length = 0
        neg_length = 0
        for result_pair in results_dict.keys():
            if args.threshold == None or args.threshold <= results_dict[result_pair]:
                if result_pair in ppis:
                    #tmp_outfilepos.write(result_pair + '\t' + str(results_dict[result_pair]) +"\n")


                    true_array.append(1)
                    prob_array.append(results_dict[result_pair])
                    true_length = true_length + 1
                elif result_pair in neg_ppis:
                    #tmp_outfileneg.write(result_pair + '\t' + str(results_dict[result_pair]) +"\n")
                    true_array.append(-1)
                    prob_array.append(results_dict[result_pair])
                    neg_length = neg_length + 1

        print "Num Trues"
        print str(true_length)
        print "Num Negs"
        print str(neg_length)
        print len(true_array)
        print len(prob_array)

        #tmp_outfileneg.close()
        #tmp_outfileneg.close()
        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        #average_precision = average_precision_score(true_array, prob_array)

        #print "precision: %s" % (list(precision),)
        #print "recall: %s" % (list(recall),)
        print len(precision)
        print len(recall)
        print len(thresholds)


        if args.plot_thresholds:
		indexOf1 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.1) 
		indexOf2 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.2) 
		indexOf3 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.3) 
		indexOf4 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.4) 
		indexOf5 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.5) 
		indexOf6 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.6) 
		indexOf7 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.7) 
		indexOf8 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.8) 
		indexOf9 = next(x[0] for x in enumerate(thresholds) if x[1] >= 0.9) 
		#indexOf10 = next(x[0] for x in enumerate(thresholds) if x[1] >= 1.0) 

		print thresholds[indexOf1]
		print indexOf1
		print thresholds[indexOf2]
		print indexOf2
		print thresholds[indexOf3]
		print indexOf3
		print thresholds[indexOf4]
		print indexOf4
		#print thresholds[indexOf10]
		#print indexOf10

		#print thresholds[indexOf1-1]
		#print thresholds[indexOf10-1]

		threshold_precisions = []
		threshold_recalls = []
		print "threshold 0.1"
		print "precision: %s" % precision[indexOf1]
		print "recall: %s" % recall[indexOf1]
		threshold_precisions.append(precision[indexOf1])
		threshold_recalls.append(recall[indexOf1])
		
		print "threshold 0.2"
		print "precision: %s" % precision[indexOf2]
		print "recall: %s" % recall[indexOf2]
		threshold_precisions.append(precision[indexOf2])
		threshold_recalls.append(recall[indexOf2])
		
		print "threshold 0.3"
		print "precision: %s" % precision[indexOf3]
		print "recall: %s" % recall[indexOf3]
		threshold_precisions.append(precision[indexOf3])
		threshold_recalls.append(recall[indexOf3])
		
		print "threshold 0.4"
		print "precision: %s" % precision[indexOf4]
		print "recall: %s" % recall[indexOf4]
		threshold_precisions.append(precision[indexOf4])
		threshold_recalls.append(recall[indexOf4])
		
		print "threshold 0.5"
		print "precision: %s" % precision[indexOf5]
		print "recall: %s" % recall[indexOf5]
		threshold_precisions.append(precision[indexOf5])
		threshold_recalls.append(recall[indexOf5])
		
		print "threshold 0.6"
		print "precision: %s" % precision[indexOf6]
		print "recall: %s" % recall[indexOf6]
		threshold_precisions.append(precision[indexOf6])
		threshold_recalls.append(recall[indexOf6])
		
		print "threshold 0.7"
		print "precision: %s" % precision[indexOf7]
		print "recall: %s" % recall[indexOf7]
		threshold_precisions.append(precision[indexOf7])
		threshold_recalls.append(recall[indexOf7])

		print "threshold 0.8"
		print "precision: %s" % precision[indexOf8]
		print "recall: %s" % recall[indexOf8]
		threshold_precisions.append(precision[indexOf8])
		threshold_recalls.append(recall[indexOf8])
		
		print "threshold 0.9"
		print "precision: %s" % precision[indexOf9]
		print "recall: %s" % recall[indexOf9]
		threshold_precisions.append(precision[indexOf9])
		threshold_recalls.append(recall[indexOf9])

		#print "threshold 1.0"
		#print "precision: %s" % precision[indexOf10]
		#print "recall: %s" % recall[indexOf10]
		#threshold_precisions.append(precision[indexOf10])
		#threshold_recalls.append(recall[indexOf10])


        label = args.results_wprob[i]
        if args.labels != None:
            label = args.labels[i]
        line, = plt.plot(recall, precision, label=label)
        if args.plot_thresholds:
            plt.scatter(threshold_recalls, threshold_precisions, color=line.get_color())


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


