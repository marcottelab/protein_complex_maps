
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
    parser.add_argument("--output_file", action="store", dest="output_file", required=False, default=None,
                                    help="Filename of output clusters for best set of parameters")
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
    parser.add_argument("--clusterone_max_overlap", action="store", dest="clusterone_max_overlap", nargs='+', required=False, 
                                    default=[0.5,0.75,0.9], 
                                    help="ClusterOne max-overlap parameter sweep, default = 0.5 0.75 0.9")
    parser.add_argument("--clusterone_seed_method", action="store", dest="clusterone_seed_method", nargs='+', required=False, 
                                    default=['nodes'], 
                                    help="ClusterOne seed method parameter sweep (nodes, cliques, unused_nodes, edges, default = nodes")
    parser.add_argument("--ppi_fraction", action="store", dest="ppi_fraction", nargs='+', required=False, 
                                    default=[0.005,0.01,0.05,.1,.25,.5,.75,1.0], 
                                    help="Use top fraction for further clustering, default = 0.005 0.01 0.05 .1 .25 .5 .75 1.0")
    parser.add_argument("--mcl", action="store", dest="mcl_bin", required=False, 
                                    default='mcl',
                                    help="Location of mcl binary, default = 'mcl' ")
    parser.add_argument("--mcl_inflation", action="store", dest="mcl_inflation", nargs='+', required=False, 
                                    default=[None],
                                    help="MCL Inflation (-I) parameter, default = [None] (no 2-stage clustering), docs suggest = 1.2 - 5.0")
    parser.add_argument("--eval_metric", action="store", dest="eval_metric", required=False, default='mmr',
                                    help="Evaluation metric used to determine best set of parameters (mmr, acc, sensitivity, ppv), default=mmr")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    parser.add_argument("--temp_dir", action="store", dest="temp_dir", required=False, default=None,
                                    help="Where to store temporary files generated during processing, default=None (defaults to OS level tmp), in memory suggestion = /dev/shm/")

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
        input_network_list = input_network_file.readlines()


    ppi_scores = dict()
    for ppi in input_network_list:
        ppi_scores[frozenset([ppi.split()[0],ppi.split()[1]])] = float(ppi.split()[2])

    random.seed(args.random_seed)
    #for i in range(args.bootstrap_iter):

    #kdrew: size and density are clusterOne parameters
    #size_sweep = [2,]
    #density_sweep=[.1,.25,.3,.35,]
    #kdrew: fraction is the fraction of top ppis to include
    #fraction_sweep=[0.005,0.01,0.05,.1,.25,.5,.75,1.0]

    size_sweep = args.clusterone_size
    overlap_sweep = args.clusterone_max_overlap
    seed_method_sweep = args.clusterone_seed_method
    density_sweep = args.clusterone_density
    fraction_sweep = args.ppi_fraction
    inflation_sweep = args.mcl_inflation

    best_size = None
    best_density = None
    best_fraction = None
    best_overlap = None
    best_seed_method = None
    best_inflation = None
    best_eval = None

    p = mp.Pool(args.procs)
    network_input_list = []
    for ii, parameters  in enumerate(it.product(size_sweep, density_sweep, fraction_sweep, overlap_sweep, seed_method_sweep, inflation_sweep )):
        #print parameters
        #kdrew: unpack parameters
        size, density, fraction, overlap, seed_method, inflation  = parameters

        #kdrew: only take the topN ppis, assumes already sorted network (probably stupidly)
        sizeOfTopNetwork = int(len(input_network_list)*float(fraction))
        network_list = input_network_list[0:sizeOfTopNetwork]

        #kdrew: cluster network_list
        parameter_dict = dict()
        parameter_dict['network_list'] = network_list
        parameter_dict['ppi_scores'] = ppi_scores
        parameter_dict['args'] = args
        parameter_dict['size'] = str(size)
        parameter_dict['density'] = str(density)
        parameter_dict['overlap'] = str(overlap)
        parameter_dict['seed_method'] = str(seed_method)
        parameter_dict['fraction'] = str(fraction)
        parameter_dict['inflation'] = str(inflation)
        parameter_dict['i'] = ii

        #kdrew: append to list of parameters for input into multiprocessor map
        network_input_list.append(parameter_dict)

    #kdrew: call clustering on parameter combinations
    cluster_predictions = p.map(cluster_helper, network_input_list)

    #kdrew: iterate through clustering results, returns tuple of prediction and number of iteration (ii)
    for cluster_prediction, ii in cluster_predictions:

        network_list = network_input_list[ii]['network_list']
        size = network_input_list[ii]['size']
        density = network_input_list[ii]['density']
        overlap = network_input_list[ii]['overlap']
        seed_method = network_input_list[ii]['seed_method']
        fraction = network_input_list[ii]['fraction']
        inflation = network_input_list[ii]['inflation']

        #kdrew: compare gold standard vs predicted clusters
        cplx_comparison = cc.ComplexComparison(gold_standard_complexes, cluster_prediction) 

        metric_dict = dict()
        metric_dict['acc'] = cplx_comparison.acc() 
        metric_dict['sensitivity'] = cplx_comparison.sensitivity() 
        metric_dict['ppv'] = cplx_comparison.ppv() 
        metric_dict['mmr'] = cplx_comparison.mmr() 
        print "size %s, density %s, overlap %s, seed_method %s, fraction %s, inflation %s, acc %s, sensitivity %s, ppv %s, mmr %s" % (size, density, overlap, seed_method, fraction, inflation, metric_dict['acc'], metric_dict['sensitivity'], metric_dict['ppv'], metric_dict['mmr'])



        #kdrew: bootstrap by resampling portions of network
        sizeOfBootstrap = int(len(network_list)*args.bootstrap_fraction)
        #bootstrapped_networks = [random.sample(network_list,sizeOfBootstrap) for i in range(args.bootstrap_iter)]
        #kdrew: generate k shuffled networks
        shuffled_networks = [ sorted( network_list, key=lambda k: random.random() ) for i in range(args.bootstrap_iter)]
        #kdrew: split shuffled networks into partitions based on size of bootstrap
        bootstrapped_networks = [ sn[:sizeOfBootstrap] for sn in shuffled_networks ]
        bootstrapped_test_networks = [ sn[sizeOfBootstrap:] for sn in shuffled_networks ]

        #multiproc_input = [(boot_net, args, str(size), str(density), str(overlap), str(seed_method), str(fraction), i) for i, boot_net in enumerate(bootstrapped_networks)]
        multiproc_input = []

        for i, boot_net in enumerate(bootstrapped_networks):
            parameter_dict = dict()
            parameter_dict['ppi_scores'] = ppi_scores
            parameter_dict['network_list'] = boot_net
            parameter_dict['args'] = args
            parameter_dict['size'] = str(size)
            parameter_dict['density'] = str(density)
            parameter_dict['overlap'] = str(overlap)
            parameter_dict['seed_method'] = str(seed_method)
            parameter_dict['fraction'] = str(fraction)
            parameter_dict['inflation'] = str(inflation)
            parameter_dict['i'] = i
            multiproc_input.append(parameter_dict)

        bootstrapped_cluster_predictions = p.map(cluster_helper, multiproc_input)

        #kdrew: compare full clustered set vs bootstrapped clusters
        multiproc_input = [(cluster_prediction, predicted_clusters, bootstrapped_test_networks[i]) for predicted_clusters, i in bootstrapped_cluster_predictions]
        bootstrap_cplx_cmp_metrics = p.map(comparison_helper, multiproc_input) 
        for boot_cmp in bootstrap_cplx_cmp_metrics:
            print "bootstrapped: size %s, density %s, overlap %s, seed_method %s, fraction %s, inflation %s,  acc %s, sensitivity %s, ppv %s, mmr %s, ppi_recovered %s" % (size, density, overlap, seed_method, fraction, inflation, boot_cmp['acc'], boot_cmp['sensitivity'], boot_cmp['ppv'], boot_cmp['mmr'], boot_cmp['percent_ppi_recovered'])


        #kdrew: keeping track of the best parameter set
        if best_eval == None or best_eval < metric_dict[args.eval_metric]: 
            best_eval = metric_dict[args.eval_metric]
            best_size = size
            best_density = density
            best_overlap = overlap
            best_seed_method = seed_method
            best_fraction = fraction
            best_inflation = inflation
            best_cluster_prediction = cluster_prediction
            print "best size: %s density: %s overlap: %s seed_method: %s fraction: %s inflation: %s numOfClusters: %s" % (best_size, best_density, best_overlap, best_seed_method, best_fraction, best_inflation, len(best_cluster_prediction))




    ##kdrew: only take the topN ppis, assumes already sorted network (probably stupidly)
    #sizeOfTopNetwork_best = int(len(input_network_list)*float(best_fraction))
    #network_list_best = input_network_list[0:sizeOfTopNetwork_best]
    #
    ##kdrew: cluster network_list
    #multiproc_input_best = (network_list_best, args, str(best_size), str(best_density), str(best_overlap), str(best_seed_method), 0)
    #cluster_prediction_best, ii = cluster_helper( multiproc_input_best )

    
    #kdrew: output best cluster prediction
    if args.output_file != None:
        with open (args.output_file, "w") as output_file:
            for cluster in best_cluster_prediction:
                output_file.write(' '.join(cluster))
                output_file.write("\n")

        
    p.close()
    p.join()



#kdrew: helper function that compares two sets of complexes/clusters, for use in multiprocessor
#kdrew: first parameter in tuple is the gold standard and second is the predicted clusters
def comparison_helper(parameter_tuple):

    gold_std = parameter_tuple[0]
    pred_clst = parameter_tuple[1]
    #kdrew: test_net is the other half of the bootstrap which is used for testing
    test_net = parameter_tuple[2]
    ppi_recovered_count = 0
    for ppi in test_net:
        prot1 = ppi.split()[0]
        prot2 = ppi.split()[1]
        for clst in pred_clst:
            if prot1 in clst and prot2 in clst:
                ppi_recovered_count += 1
                #kdrew: only count the ppi once if found in cluster
                break

    cplx_cmp = cc.ComplexComparison(gold_std, pred_clst) 
    d = dict()
    d['acc'] = cplx_cmp.acc()
    d['sensitivity'] = cplx_cmp.sensitivity()
    d['ppv'] = cplx_cmp.ppv()
    d['mmr'] = cplx_cmp.mmr()
    d['percent_ppi_recovered'] = (1.0*ppi_recovered_count) / len(test_net)

    return d

def cluster_helper(parameter_dict):

    input_network_list = parameter_dict['network_list']
    ppi_scores = parameter_dict['ppi_scores']
    args = parameter_dict['args']
    size = parameter_dict['size']
    density = parameter_dict['density']
    overlap = parameter_dict['overlap']
    seed_method = parameter_dict['seed_method']
    fraction = parameter_dict['fraction']
    i = parameter_dict['i']
    inflation = parameter_dict['inflation']

    #kdrew: create temp file for bootstrapped input network, clusterone requires a file input
    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
    try:
        #kdrew: write input_network or bootstrapped input network to temp file
        #fileTemp.write(input_network_str)
        for bstrap_ppi in input_network_list:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        #kdrew: run clustering
        proc = sp.Popen(['java', '-jar', args.clustone_jar, fileTemp.name, '-s', size, '-d', density, '--max-overlap', overlap, '--seed-method', seed_method], stdout=sp.PIPE, stderr=sp.PIPE)
        clustone_out, err = proc.communicate()

    finally:
        os.remove(fileTemp.name)

    #kdrew: take output of cluster one (predicted clusters) and store them into list
    predicted_clusters = []
    for line in clustone_out.split('\n'):
        if len(line.split() ) > 0:
            predicted_clusters.append(line.split())

    #kdrew: converted inflation into str for use on commandline but if None check against str(None)
    if inflation != str(None):
        mcl_clusters = []
        #kdrew: for each predicted cluster, recluster using MCL
        for clust in predicted_clusters:
            #kdrew: for every pair in cluster find edge weight in input_network_list(?)
            fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
            outTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
            #print fileTemp.name
            #print outTemp.name
            try:
                for prot1, prot2 in it.combinations(clust,2):
                    try:
                        ppi_score = ppi_scores[frozenset([prot1,prot2])] 
                    except KeyError:
                        ppi_score = 0.0

                    ppi_str = "%s\t%s\t%s\n" % (prot1, prot2, ppi_score)
                    fileTemp.write(ppi_str)
                fileTemp.close()

                proc = sp.Popen([args.mcl_bin, fileTemp.name, '--abc', '-o', outTemp.name, '-I', inflation], stdout=sp.PIPE, stderr=sp.PIPE)
                mcl_out, err = proc.communicate()

                #print mcl_out
                #print err

                outfile = open(outTemp.name,"rb")
                for line in outfile.readlines():
                    mcl_clusters.append(line.split())
                outfile.close()


            finally:
                os.remove(fileTemp.name)
                os.remove(outTemp.name)
                #print "in finally"

        predicted_clusters = mcl_clusters


    return predicted_clusters, i


if __name__ == "__main__":
    main()


