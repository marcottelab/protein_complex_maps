
import numpy as np
#kdrew: doesn't need matplotlib, remove (?)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import re
import operator

import argparse
import pickle
import pandas as pd
import itertools as it

import seaborn as sns

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
    parser.add_argument("--plot_heatmap", action="store_true", dest="plot_heatmap", required=False, default=False,
                                    help="Plot heatmap of c and g values")
    parser.add_argument("--threshold_recall", action="store", type=float, dest="threshold_recall", required=False, 
                                    help="Only calculate predictions at recall threshold")
    parser.add_argument("--parse_filename", action="store_true", dest="parse_filename", required=False, default=False,
            help="Search filename for parameters and leaveout set id, (example results: leaveout0 and c2.0_g0.5, example ppis: leaveout_ppis0), default=False")

    args = parser.parse_args()

    outfile = open(args.output_file,"wb")

    leaveout_regex = r"leaveout[0-9]+"
    leaveout_ppis_regex = r"leaveout_ppis[0-9]+"
    c_regex = r"c[0-9]+\.[0-9]+"
    g_regex = r"g[0-9]+\.[0-9]+"
    cg_dict = dict()
    pos_dict = dict()
    neg_dict = dict()
    if args.parse_filename:
        for rfile in args.results_wprob:
            leaveout_match = re.findall(leaveout_regex, rfile)[0]
            leaveout_id = leaveout_match[-1]
            c_match = re.findall(c_regex, rfile)[0]
            g_match = re.findall(g_regex, rfile)[0]
            try:
                cg_dict[(c_match, g_match)][leaveout_id] = rfile
            except KeyError:
                cg_dict[(c_match, g_match)] = dict()
                cg_dict[(c_match, g_match)][leaveout_id] = rfile

        for pfile in args.positives:
            leaveout_ppis_match = re.findall(leaveout_ppis_regex, pfile)[0]
            leaveout_id = leaveout_ppis_match[-1]
            pos_dict[leaveout_id] = pfile
        for nfile in args.negatives:
            leaveout_neg_ppis_match = re.findall(leaveout_ppis_regex, nfile)[0]
            leaveout_id = leaveout_neg_ppis_match[-1]
            neg_dict[leaveout_id] = nfile

        #print cg_dict
        #print pos_dict
        #print neg_dict

        cg_results = dict()
        for cg_pair in cg_dict:
            results_wprob_filenames = []
            positives_filenames = []
            negatives_filenames = []
            for k in cg_dict[cg_pair]:
                results_wprob_filenames.append(cg_dict[cg_pair][k])
                positives_filenames.append(pos_dict[k])
                negatives_filenames.append(neg_dict[k])


            cg_results[cg_pair] = calc_metric_files(results_wprob_filenames, positives_filenames, negatives_filenames)

        for item in sorted(cg_results.items(), key=operator.itemgetter(1)):
            k = item[0]
            outfile.write("cg: %s mean pr auc: %s\n" % (k, cg_results[k]) )
    else:
        outfile.write("mean pr auc: %s\n" % calc_metric_files(args.results_wprob, args.positives, args.negatives))

    if args.plot_heatmap:
        plot_heatmap(cg_results, args.output_file)

def plot_heatmap(cg_results, output_file):
    c_vals = set()
    g_vals = set()
    for cg in cg_results.keys():
        #kdrew: ugly parsing to remove 'c' or 'g' of keys and turn into floats
        c_vals.add(float(cg[0][1:]))
        g_vals.add(float(cg[1][1:]))

    c_vals = sorted(c_vals)
    g_vals = sorted(g_vals)

    results = []
    for c in c_vals:
        results2 = []
        for g in g_vals:
            try:
                #kdrew: add 'c' and 'g' back into values for referencing results in dictionary
                results2.append(cg_results[("c%s"%c,"g%s"%g)])
            except KeyError:
                results2.append(np.nan)
        results.append(results2)
    ax = sns.heatmap(results, xticklabels=g_vals, yticklabels=c_vals, cbar_kws={'label': 'Precision Recall AUC'})
    ax.set(xlabel='gamma', ylabel='C')
    plt.tight_layout()
    
    ax.get_figure().savefig("%s.pdf" % output_file)

def calc_metric_files(results_wprob_filenames, positives_filenames, negatives_filenames):

    results_dict_list = []
    all_proteins_list = []
    ppi_list = []
    neg_ppi_list = []
    for i, filename in enumerate(results_wprob_filenames):
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


        positive_file = open(positives_filenames[i], "rb")
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
        if negatives_filenames == None:
            for cpair in it.combinations(all_proteins,2):
                pair = str(sorted(list(frozenset(cpair))))
                if pair not in ppis:
                    #print "pair is neg: %s" % ' '.join(pair)
                    neg_ppis.add(pair)
        else:
            #kdrew: readin from file
            negative_file = open(negatives_filenames[i],"rb")
            for line in negative_file.readlines():
                if len(line.split()) >= 2:
                    id1 = line.split()[0]
                    id2 = line.split()[1]
                    ppi_str = str(sorted(list(frozenset([id1,id2]))))
                    neg_ppis.add(ppi_str)

        neg_ppi_list.append(neg_ppis)

    return calc_mean_pr_auc(results_dict_list, ppi_list, neg_ppi_list)

def calc_mean_pr_auc(results_dict_list, ppi_list, neg_ppi_list):
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

        precision, recall, thresholds = precision_recall_curve(true_array, prob_array) 
        aps = average_precision_score(true_array, prob_array)
        prauc_list.append(aps)

        #print "pr_auc: %s" % (aps,)

    return np.mean(prauc_list)


if __name__ == "__main__":
    main()


