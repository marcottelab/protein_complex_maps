
import os
import argparse
import pickle
import random
import numpy as np
import pandas as pd
import itertools as it
import subprocess as sp
import multiprocessing as mp
import tempfile as tf

import protein_complex_maps.complex_comparison as cc

def main():

    parser = argparse.ArgumentParser(description="Finds optimal parameters for clustering of ppi network")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True, 
                                    help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True, 
                                    help="Filename of gold standard complexes, one cluster per line, ids tab-separated")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", required=False, default=None,
                                    help="Sets random seed (int), default=None")
    parser.add_argument("--bootstrap_iter", action="store", type=int, dest="bootstrap_iter", required=False, default=10,
                                    help="Number of bootstrap iterations (int, default=10)")
    parser.add_argument("--bootstrap_fraction", action="store", type=float, dest="bootstrap_fraction", required=False, default=0.5,
                                    help="Fraction of edges to sample for bootstrapping, default=0.5")
    parser.add_argument("--clusterone", action="store", dest="clustone_jar", required=False, 
                                    default="/home/kdrew/programs/clusterone/cluster_one-1.0.jar",
                                    help="Location of cluster_one jar file, default= /home/kdrew/programs/clusterone/cluster_one-1.0.jar")
    parser.add_argument("--clusterone_size", action="store", dest="clusterone_size", nargs='+', required=False, 
                                    default=[2,],
                                    help="ClusterOne Size parameter sweep, default = 2")
    parser.add_argument("--clusterone_density", action="store", dest="clusterone_density", nargs='+', required=False, 
                                    default=[.1,.25,.3,.35,], 
                                    help="ClusterOne Density parameter sweep, default = .1 .25 .3 .35")
    parser.add_argument("--ppi_fraction", action="store", dest="ppi_fraction", nargs='+', required=False, 
                                    default=[0.005,0.01,0.05,.1,.25,.5,.75,1.0], 
                                    help="Use top fraction for further clustering, default = 0.005 0.01 0.05 .1 .25 .5 .75 1.0")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")

    args = parser.parse_args()

    #kdrew: read gold standard into list
    gstd_file = open(args.gold_standard_filename,"rb")
    gold_standard_complexes = []
    for line in gstd_file.readlines():
        gold_standard_complexes.append(line.split())


    #kdrew: original commandline for clusterone
    #java -jar ~/programs/clusterone/cluster_one-1.0.jar blake_bioplex_merge_wkeys_deduped_corum_train_labeled.libsvm0.scaleByTrain.resultWprob_pairs_noself_nodups_wprob.txt > blake_bioplex_merge_wkeys_deduped_corum_train_labeled.libsvm0.scaleByTrain.resultWprob_pairs_noself_nodups_wprob.txt.clusterone

    #kdrew: read in the input network into a string
    with open (args.input_network, "r") as input_network_file:
        #data=input_network_file.read().replace('\n', '')
        #input_network_str = input_network_file.read()
        input_network_list = input_network_file.readlines()

    random.seed(args.random_seed)
    #for i in range(args.bootstrap_iter):

    #kdrew: size and density are clusterOne parameters
    #size_sweep = [2,]
    #density_sweep=[.1,.25,.3,.35,]
    #kdrew: fraction is the fraction of top ppis to include
    #fraction_sweep=[0.005,0.01,0.05,.1,.25,.5,.75,1.0]

    size_sweep = args.clusterone_size
    density_sweep = args.clusterone_density
    fraction_sweep = args.ppi_fraction


    p = mp.Pool(args.procs)
    for size, density, fraction in it.product(size_sweep, density_sweep, fraction_sweep):
        multiproc_input = [(gold_standard_complexes, input_network_list, args, str(size), str(density), str(fraction))]*args.bootstrap_iter
        bootstrapping_results = p.map(bootstrap_helper, multiproc_input)
        mean_acc = np.mean([result.acc() for result in bootstrapping_results])
        mean_sensitivity = np.mean([result.sensitivity() for result in bootstrapping_results])
        mean_ppv = np.mean([result.ppv() for result in bootstrapping_results])
        mean_mmr = np.mean([result.mmr() for result in bootstrapping_results])
        print "size %s, density %s, fraction %s, mean_acc %s, mean_sensitivity %s, mean_ppv %s, mean_mmr %s" % (size, density, fraction, mean_acc, mean_sensitivity, mean_ppv, mean_mmr)

    p.close()
    p.join()

def bootstrap_helper(parameter_tuple):
    gold_standard_complexes = parameter_tuple[0]
    input_network_list = parameter_tuple[1]
    args = parameter_tuple[2]
    size = parameter_tuple[3]
    density = parameter_tuple[4]
    fraction = float(parameter_tuple[5])

    #kdrew: only take the topN ppis, assumes already sorted network (probably stupidly)
    sizeOfTopNetwork = int(len(input_network_list)*fraction)
    network_list = input_network_list[0:sizeOfTopNetwork]

    #kdrew: bootstrap by resampling portions of network
    bootstrap_fraction = args.bootstrap_fraction
    sizeOfBootstrap = int(len(network_list)*bootstrap_fraction)
    bootstrapped_network = random.sample(network_list,sizeOfBootstrap)
    #kdrew: create temp file for bootstrapped input network, clusterone requires a file input
    fileTemp = tf.NamedTemporaryFile(delete=False)
    try:
        #kdrew: write input_network or bootstrapped input network to temp file
        #fileTemp.write(input_network_str)
        for bstrap_ppi in bootstrapped_network:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        #kdrew: run clustering
        proc = sp.Popen(['java', '-jar', args.clustone_jar, fileTemp.name, '-s', size, '-d', density], stdout=sp.PIPE, stderr=sp.PIPE)
        clustone_out, err = proc.communicate()

    finally:
        os.remove(fileTemp.name)

    #kdrew: take output of cluster one (predicted clusters) and store them into list
    predicted_clusters = []
    for line in clustone_out.split('\n'):
        if len(line.split() ) > 0:
            predicted_clusters.append(line.split())


    #kdrew: compare gold standard vs predicted clusters
    cplx_compare = cc.ComplexComparison(gold_standard_complexes, predicted_clusters)
    #print "Sensitivity: %s" % cplx_compare.sensitivity()
    #print "PPV: %s" % cplx_compare.ppv()
    #print "ACC: %s" % cplx_compare.acc()

    return cplx_compare


if __name__ == "__main__":
    main()


