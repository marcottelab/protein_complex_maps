    
import numpy as np
import numpy.random as rand
import pandas as pd
import pickle as p
import argparse
import itertools as it
import random
import bisect
import scipy.misc as misc
from scipy.stats import hmean


import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import protein_complex_maps.complex_comparison as cc
import protein_complex_maps.complexes.noise_complexes as nc



def main():

    parser = argparse.ArgumentParser(description="Compare cluster predictions to gold standard complexes")
    parser.add_argument("--cluster_predictions", action="store", dest="cluster_filename", required=True, 
                                            help="Filename of cluster predictions, format one cluster per line, ids space separated")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True, 
                                            help="Filename of gold standard complexes, format one complex per line, ids space separated")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename for plotting histograms of metrics")
    parser.add_argument("--shuffle_fraction", action="store", dest="shuffle_fraction", type=float, required=False, default=0.0,
                                            help="Fraction of complex memberships to shuffle, range: 0.0 - 1.0")
    parser.add_argument("--iterations", action="store", dest="iterations", type=int, required=False, default=10,
                                            help="Number of random samples to generate")

    args = parser.parse_args()

    gold_standard_complexes = []
    gold_file = open(args.gold_standard_filename,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())
    gold_file.close()

    predicted_clusters = []
    clpred_f = open(args.cluster_filename,"rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())
    clpred_f.close()


    for i in xrange(args.iterations):
        #noise_complexes = nc.shuffle_complexes(predicted_clusters, args.shuffle_fraction)
        #noise_complexes = nc.breakup_complexes(predicted_clusters, args.shuffle_fraction)
        #noise_complexes = nc.remove_by_size_complexes(predicted_clusters, upper_size_threshold=5)
        #noise_complexes = nc.remove_by_size_complexes(predicted_clusters, lower_size_threshold=5)
        #noise_complexes = nc.powerset_complexes(predicted_clusters) 
        #noise_complexes = nc.add_random_complexes(predicted_clusters) 
        noise_complexes = nc.remove_fraction_complexes(predicted_clusters, remove_upper=False) 

        cplx_compare = cc.ComplexComparison(gold_standard_complexes, noise_complexes)
        print "Sensitivity: %s" % cplx_compare.sensitivity()
        print "PPV: %s" % cplx_compare.ppv()
        print "ACC: %s" % cplx_compare.acc()
        print "MMR: %s" % cplx_compare.mmr()
        print "PWMMR: %s" % cplx_compare.pwmmr()
        print "MMR_PWMMR_hmean: %s" % cplx_compare.mmr_pwmmr_hmean()
        ccmm = cplx_compare.clique_comparison_metric_mean()
        print "Clique Precision Mean: %s Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
        ccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
        print "Clique Weighted Precision Mean: %s Weighted Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
        print "Clique Weighted hmean: %s" % (hmean([ccmm['precision_mean'],ccmm['recall_mean']]))


    if args.plot_filename != None:
        f, subplots = plt.subplots(3)
        subplots[0].hist(cplx_compare.max_matching_ratio_distribution())
        subplots[0].set_title('MMR')
        subplots[1].hist(cplx_compare.sensitivity_distribution())
        subplots[1].set_title('Sensitivity')
        subplots[2].hist(cplx_compare.ppv_distribution())
        subplots[2].set_title('PPV')
        plt.savefig(args.plot_filename)

if __name__ == "__main__":
        main()

