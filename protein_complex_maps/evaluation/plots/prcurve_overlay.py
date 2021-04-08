
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import pickle
import pandas as pd
import itertools as it

from collections import Counter #for averaging dictionary values
import pandas as pd #for averaging dictionary values

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import average_precision_score


from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
#rcParams.update({'grid': False}) # KeyError: u'grid is not a valid rc parameter.See rcParams.keys() for a list of valid parameters.'
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

def main():

    parser = argparse.ArgumentParser(description="Generate Precision Recall Curve")
    parser.add_argument("--results_wprob", action="store", dest="results_wprob", nargs='+', required=True, 
                                    help="Filename of pair results with probability")
    parser.add_argument("--input_positives", action="store", dest="positives", nargs='+', required=True, 
                                    help="Filename of positive pairs")
    parser.add_argument("--input_negatives", action="store", dest="negatives", nargs='+',
                                    help="Filename of negative pairs, default = None (generated from processing positives)")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--labels", action="store", dest="labels", nargs='+', required=False, default=None,
                                    help="Labels for input results in order of sets of pairs")
    parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, 
                                    help="Only tally predictions above probability threshold")
    parser.add_argument("--plot_thresholds", action="store_true", dest="plot_thresholds", required=False, default=False,
                                    help="Add probability threshold markers to plot")
    parser.add_argument("--average_probs", action="store_true", dest="avg_probs", required=False, default=False,
                                    help="Average a set of result_wprob instead of plotting individually")

    parser.add_argument("--header", action="store_true", dest="header", required=False, default=False,
                                    help="Set true in results_wprob has a header")

    args = parser.parse_args()

    output_pr_table = pd.DataFrame()
    output_roc_table = pd.DataFrame()


    results_dict_list = []
    for filename in args.results_wprob:
        print(filename)
        results_dict = dict()
        results_file = open(filename,"rb")
        for line in results_file.readlines():
            try: 
                if len(line.split()) == 3:
                    id1 = line.split()[0]
                    id2 = line.split()[1]
                    ppi_str = str(sorted(list(frozenset([id1,id2]))))
                    results_dict[ppi_str] = float(line.split()[2])
            except:
                print "WARNING: did not read line: %s" % (line,)

        print "results size: %s" % (len(results_dict))
        results_dict_list.append(results_dict)

    if args.avg_probs:
        df = pd.DataFrame(results_dict_list)
        print(df)
        results_dict = dict(df.mean())
        print(results_dict)
        results_dict_list = []
        results_dict_list.append(results_dict)
         
    #assert len(args.positives) == len(args.labels)

   
    for evalset in range(0, len(args.positives)):
        all_proteins = set()
        ppis = set()
        positive_file = open(args.positives[evalset],"rb")
        print(args.positives[evalset])
        for line in positive_file.readlines():
                if len(line.split()) >= 2:
                    id1 = line.split()[0]
                    id2 = line.split()[1]
                    all_proteins.add(id1)
                    all_proteins.add(id2)
        
                    ppi_str = str(sorted(list(frozenset([id1,id2]))))
                    ppis.add(ppi_str)
 

        for evalset2 in range(0,len(args.negatives)):

            negative_file = open(args.negatives[evalset2],"rb")
            neg_ppis = set()

            for line in negative_file.readlines():
                    if len(line.split()) >= 2:
                        id1 = line.split()[0]
                        id2 = line.split()[1]
                        ppi_str = str(sorted(list(frozenset([id1,id2]))))
                        neg_ppis.add(ppi_str)
 

  
       
        
            print(args.negatives[evalset2])        
            for i, results_dict in enumerate(results_dict_list):
                print(i)
                try:
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

                    #print "FoundNegs: %s" % (len([i for i in true_array if i == -1]))
                    #print "FoundPos: %s" % (len([i for i in true_array if i == 1]))


                    print "#ofNegs: %s" % (len(sorted_neg_list))
                    print "#ofPos: %s" % (len(sorted_pos_list))
                    print "FDR: %s" % (1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))
            
                    print "#ofNegs: %s, #ofPos: %s, FDR: %s" % (len(sorted_neg_list),len(sorted_pos_list),1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))
                    #for k in sorted_neg_list[-100:]:
                    #    print "neg: %s:%s" % (k[0], k[1])
                    #for k in sorted_pos_list[-100:]:
                    #    print "pos: %s:%s" % (k[0], k[1])
                    print len(true_array)
                    print len(prob_array)
                                
                    precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
                    print(thresholds)
                    print("pr dims")
                    print(len(recall))
                    print(len(precision))
                    print(len(thresholds))
                    print(recall)
                    print(precision)
                    print('pr thresholds', thresholds)
                    FPR, TPR, roc_thresholds = roc_curve(true_array, prob_array)
                    print("roc dims")
                    print len(FPR)
                    print len(TPR) 
                    print len(roc_thresholds)
                    print('roc_thresholds', roc_thresholds)
                    #average_precision = average_precision_score(true_array, prob_array)
           
          
 
                    #print "precision: %s" % (list(precision),)
                    #print "recall: %s" % (list(recall),)
                    print len(precision)
                    print len(recall)
                    print len(thresholds)
            
            
                    label = args.results_wprob[i]
                    
                    if args.labels != None:
                        #label = args.labels[evalset]
                        label = args.labels[i] + "_" + str(evalset) + "_" + str(evalset2)

                    if i == 0 and evalset == 0 and evalset2 ==0:
                          label = args.labels[0]

                    elif i == 0 and evalset == 1 and evalset2 ==1:
                          label = args.labels[1]
                          print(precision)
                          print(recall)
                    #elif i == 1 and evalset == 0 and evalset2 ==0:
                    #      label = 'cross_val'
                    else:
                         print("is this happening?", i, evalset, evalset2) 
                         continue
                    print("TESTING", i, evalset, evalset2)
                    print("recall")
                    print(recall)
                    print("precision")
                    print(precision)
                    print("label")
                    print(label) 

                    print "#ofNegs: %s" % (len(sorted_neg_list))
                    print "#ofPos: %s" % (len(sorted_pos_list))
                    print "FDR: %s" % (1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))
            
                    print "#ofNegs: %s, #ofPos: %s, FDR: %s" % (len(sorted_neg_list),len(sorted_pos_list),1.0*len(sorted_neg_list)/(len(sorted_pos_list)+len(sorted_neg_list)))

                    label_pr_table = pd.DataFrame({'Recall':recall, 'Precision':precision, 'label':label})
                    
                    print(label_pr_table)
                    print("HERE")
                    output_pr_table = output_pr_table.append(label_pr_table, ignore_index=True)
                    
                    label_roc_table = pd.DataFrame({'TPR':TPR, 'FPR':FPR,'Thresholds':roc_thresholds, 'label':label})
                    print(label_roc_table)
                    output_roc_table = output_roc_table.append(label_roc_table, ignore_index=True)

                    pr_filename = label + "_pr_tresholds.csv"
                    pd.DataFrame({'pr_thresholds':thresholds}).to_csv(pr_filename, index = False)

                    roc_filename = label + "_roc_tresholds.csv"
          
                    pd.DataFrame({'roc_thresholds':roc_thresholds}).to_csv(roc_filename, index = False)
 


                    line, = plt.plot(recall, precision, label=label)
                    if args.plot_thresholds:
                        plt.scatter(threshold_recalls, threshold_precisions, color=line.get_color())
                        for i, label in enumerate(threshold_labels):
                            plt.annotate(label, xy=(threshold_recalls[i],threshold_precisions[i]))
                except Exception as E:
                 print(E) 
        
    #plt.clf()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    #plt.title('Precision-Recall example: AUC={0:0.2f}'.format(average_precision))
    plt.title('Precision-Recall')
    plt.legend(loc="upper right",fontsize=8)

    plt.savefig(args.output_file)
 
    output_csv = args.output_file + ".csv"
    output_pr_table.to_csv(output_csv, index = False)

 
    output_roc_csv = args.output_file + "_roc.csv"
    output_roc_table.to_csv(output_roc_csv, index = False)

if __name__ == "__main__":
    main()


