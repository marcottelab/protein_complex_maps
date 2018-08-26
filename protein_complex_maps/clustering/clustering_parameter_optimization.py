
import os
import sys
import argparse
import pickle
import random
import numpy as np
import pandas as pd
import itertools as it
import subprocess as sp
import multiprocessing as mp
import tempfile as tf
import shutil


#kdrew: CLAIRE WTF!?!
#sys.path.append('/project/cmcwhite/protein_complex_maps/protein_complex_maps')



import networkx as nx
import agglomcluster.agglomod as ag

import protein_complex_maps.evaluation.complex_comparison as cc

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

    parser.add_argument("--ppi_fraction", action="store", dest="ppi_fraction", nargs='+', required=False, 
                                    default=[0.005,0.01,0.05,.1,.25,.5,.75,1.0], 
                                    help="Use top fraction for further clustering, default = 0.005 0.01 0.05 .1 .25 .5 .75 1.0")
    parser.add_argument("--ppi_threshold_score", action="store", dest="ppi_threshold_score", nargs='+', required=False, 
                                    default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],
                                    help="Use ppis with score higher or equal to threshold score, default = 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1")

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

    parser.add_argument("--clixo_bin", action="store", dest="clixo_bin", required=False, 
                                    default='~/programs/mhk7-clixo_0.3-a2b23b0/clixo')
    parser.add_argument("--clixo_alpha", action="store", dest="clixo_alpha", nargs='+', required=False, 
                                    default=[None],
                                    help="Clixo Alpha parameter, default = None")
    parser.add_argument("--clixo_beta", action="store", dest="clixo_beta", nargs='+', required=False, 
                                    default=[None],
                                    help="Clixo Beta parameter, default = None")

    parser.add_argument("--mcl", action="store", dest="mcl_bin", required=False, 
                                    default='mcl',
                                    help="Location of mcl binary, default = 'mcl' ")
    parser.add_argument("--mcl_inflation", action="store", dest="mcl_inflation", nargs='+', required=False, 
                                    default=[None],
                                    help="MCL Inflation (-I) parameter, default = [None] (no 2-stage clustering), docs suggest = 1.2 - 5.0")

    parser.add_argument("--cfinder", action="store", dest="cfinder_exe", required=False, 
                                    default="/home/kdrew/programs/CFinder-2.0.6--1448/CFinder_commandline64",
                                    help="Location of CFinder executable, default= /home/kdrew/programs/CFinder-2.0.6--1448/CFinder_commandline64")
    parser.add_argument("--cfinder_license", action="store", dest="cfinder_license", required=False, 
                                    default="/home/kdrew/programs/CFinder-2.0.6--1448/licence.txt",
                                    help="Location of CFinder license, default= /home/kdrew/programs/CFinder-2.0.6--1448/licence.txt")
    parser.add_argument("--cfinder_cliquesize", action="store", dest="cfinder_cliquesize", nargs='+', required=False, 
                                    default=[None],
                                    help="Cfinder clique size (-k) parameter, default = None (use CFinder's default setting, recommended: 3)")
    parser.add_argument("--cfinder_timeout", action="store", dest="cfinder_timeout", nargs='+', required=False, 
                                    default=[None],
                                    help="Cfinder timeout (-t) parameter, default = None (use CFinder's default setting, recommended: 10)")

    parser.add_argument("--trim2threshold", action="store_true", dest="trim2threshold", required=False, default=False,
                                    help="Trim final clusters of subunits that do not have an edge that passes the threshold_score, default=False")
    parser.add_argument("--twostep_combination", action="store", dest="twostep_combination", nargs='+', required=False, 
                                    default=['clusterone','mcl'],
                                    help="Combination of two step clustering, default = [clusterone,mcl], options=[clusterone,mcl,cfinder,agglomod,clixo]")

    parser.add_argument("--bootstrap_iter", action="store", type=int, dest="bootstrap_iter", required=False, default=10,
                                    help="Number of bootstrap iterations (int, default=10)")
    parser.add_argument("--bootstrap_fraction", action="store", type=float, dest="bootstrap_fraction", required=False, default=0.5,
                                    help="Fraction of edges to sample for bootstrapping, default=0.5")

    parser.add_argument("--eval_metric", action="store", dest="eval_metric", required=False, default='mmr',
                                    help="Evaluation metric used to determine best set of parameters (mmr, acc, sensitivity, ppv, clique_precision_mean, clique_recall_mean), default=mmr")
    parser.add_argument("--output_all", action="store_true", dest="output_all", required=False, default=False,
                                    help="Output all clusterings, default=False")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    parser.add_argument("--temp_dir", action="store", dest="temp_dir", required=False, default=None,
                                    help="Where to store temporary files generated during processing, default=None (defaults to OS level tmp), in memory suggestion = /dev/shm/")
    parser.add_argument("--nodelete", action="store_true", dest="nodelete", required=False, default=False,
                                    help="When set, does not delete temporary files")


    args = parser.parse_args()

    if not os.path.isfile(args.clustone_jar):
        raise IOError("File not found: %s" % args.clustone_jar)

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

    #kdrew: parse input network into dictionary of id1, id2 = score
    ppi_scores = dict()
    for ppi in input_network_list:
        ppi_scores[frozenset([ppi.split()[0],ppi.split()[1]])] = float(ppi.split()[2])

    random.seed(args.random_seed)

    size_sweep = args.clusterone_size
    overlap_sweep = args.clusterone_max_overlap
    seed_method_sweep = args.clusterone_seed_method
    density_sweep = args.clusterone_density
    fraction_sweep = args.ppi_fraction
    threshold_score_sweep = args.ppi_threshold_score
    inflation_sweep = args.mcl_inflation
    clixo_alpha_sweep = args.clixo_alpha
    clixo_beta_sweep = args.clixo_beta
    cliquesize_sweep = args.cfinder_cliquesize
    timeout_sweep = args.cfinder_timeout

    best_size = None
    best_ii = None
    best_density = None
    best_fraction = None
    best_threshold_score = None
    best_overlap = None
    best_seed_method = None
    best_inflation = None
    best_clixo_alpha = None
    best_clixo_beta = None
    best_cliquesize = None
    best_timeout = None
    best_eval = None

    p = mp.Pool(args.procs)
    network_input_list = []
    for ii, parameters  in enumerate(it.product(size_sweep, density_sweep, fraction_sweep, overlap_sweep, seed_method_sweep, inflation_sweep, clixo_alpha_sweep, clixo_beta_sweep, cliquesize_sweep, timeout_sweep, threshold_score_sweep )):
        #/print ii, parameters
        sys.stdout.flush()
        #kdrew: unpack parameters
        size, density, fraction, overlap, seed_method, inflation, clixo_alpha, clixo_beta, cliquesize, timeout, threshold_score = parameters
        threshold_score = float(threshold_score)

        if threshold_score == None:
            #print "threshold_score == None, using fraction"
            #sys.stdout.flush()

            #kdrew: only take the topN ppis, assumes already sorted network (probably stupidly)
            sizeOfTopNetwork = int(len(input_network_list)*float(fraction))
            network_list = input_network_list[0:sizeOfTopNetwork]
            #print network_list[-1]
            threshold_score = float(network_list[-1].split()[2])
            #print threshold_score
        else:
            #print "using threshold_score: %s" % (len(input_network_list))
            network_list = []
            for x in input_network_list:                   
                if float(x.split()[2]) >= float(threshold_score):
                    network_list.append(x)
                else:
                    #kdrew: assuming sorted
                    break

            fraction = 1.0*len(network_list) / len(input_network_list)

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
        parameter_dict['threshold_score'] = threshold_score
        parameter_dict['inflation'] = str(inflation)
        parameter_dict['clixo_alpha'] = str(clixo_alpha)
        parameter_dict['clixo_beta'] = str(clixo_beta)
        parameter_dict['cliquesize'] = str(cliquesize)
        parameter_dict['timeout'] = str(timeout)
        parameter_dict['twostep_combination'] = args.twostep_combination
        parameter_dict['trim2threshold'] = args.trim2threshold
        parameter_dict['i'] = ii
        parameter_dict['nodelete'] = args.nodelete 

        #kdrew: append to list of parameters for input into multiprocessor map
        network_input_list.append(parameter_dict) 
    #print network_input_list
    #kdrew: call clustering on parameter combinations
    cluster_predictions = p.map(cluster_helper, network_input_list)

    #kdrew: iterate through clustering results, returns tuple of prediction and number of iteration (ii)
    for cluster_prediction, ii in cluster_predictions:
        print "ii", ii
        network_list = network_input_list[ii]['network_list']
        size = network_input_list[ii]['size']
        density = network_input_list[ii]['density']
        overlap = network_input_list[ii]['overlap']
        seed_method = network_input_list[ii]['seed_method']
        fraction = network_input_list[ii]['fraction']
        threshold_score = network_input_list[ii]['threshold_score']
        inflation = network_input_list[ii]['inflation']
        clixo_alpha = network_input_list[ii]['clixo_alpha']
        clixo_beta = network_input_list[ii]['clixo_beta']
        timeout = network_input_list[ii]['timeout']
        cliquesize = network_input_list[ii]['cliquesize']
        twostep_combination = network_input_list[ii]['twostep_combination']
        trim2threshold = network_input_list[ii]['trim2threshold']
        nodelete = network_input_list[ii]['nodelete']

        print "ii %s, size %s, density %s, overlap %s, seed_method %s, fraction %s, threshold_score %s, inflation %s, clixo_alpha %s, clixo_beta %s, cliquesize %s, timeout %s, twostep_combination: %s, trim2threshold: %s" % (ii, size, density, overlap, seed_method, fraction, threshold_score,  inflation, clixo_alpha, clixo_beta, cliquesize, timeout, str(twostep_combination), str(trim2threshold))



        #kdrew: compare gold standard vs predicted clusters
        cplx_comparison = cc.ComplexComparison(gold_standard_complexes, cluster_prediction) 
        cplx_comparison_normalize = cc.ComplexComparison(gold_standard_complexes, cluster_prediction, normalize_by_combinations=True, pseudocount=0.00001) 

        print cplx_comparison

        try:  
            metric_dict = dict()
            metric_dict['acc'] = cplx_comparison.acc() 
            metric_dict['sensitivity'] = cplx_comparison.sensitivity() 
            metric_dict['ppv'] = cplx_comparison.ppv() 
            metric_dict['mmr'] = cplx_comparison.mmr() 
            metric_dict['precision_recall_product'] = cplx_comparison.precision_recall_product() 
            ccmm = cplx_comparison.clique_comparison_metric_mean()
            metric_dict['clique_precision_mean'] = ccmm['precision_mean']
            metric_dict['clique_recall_mean'] = ccmm['recall_mean']

            ccmm_normalize = cplx_comparison_normalize.clique_comparison_metric_mean()
            metric_dict['clique_precision_mean_normalize'] = ccmm_normalize['precision_mean']
            metric_dict['clique_recall_mean_normalize'] = ccmm_normalize['recall_mean']
            ccmm_normalize_weighted = cplx_comparison_normalize.clique_comparison_metric_mean(weighted=True)
            metric_dict['clique_precision_mean_normalize_weighted'] = ccmm_normalize_weighted['precision_mean']
            metric_dict['clique_recall_mean_normalize_weighted'] = ccmm_normalize_weighted['recall_mean']
        except Exception as e:
            print e
            continue

        print "ii %s, size %s, density %s, overlap %s, seed_method %s, fraction %s, threshold_score %s, inflation %s, clixo_alpha %s, clixo_beta %s,cliquesize %s, timeout %s, twostep_combination: %s, trim2threshold: %s, acc %s, sensitivity %s, ppv %s, mmr %s, precision_recall_product %s, clique_precision_mean %s, clique_recall_mean %s, clique_precision_mean_normalize %s, clique_recall_mean_normalize %s, clique_precision_mean_normalize_weighted %s, clique_recall_mean_normalize_weighted %s" % (ii, size, density, overlap, seed_method, fraction, threshold_score,  inflation, clixo_alpha, clixo_beta, cliquesize, timeout, str(twostep_combination), str(trim2threshold), metric_dict['acc'], metric_dict['sensitivity'], metric_dict['ppv'], metric_dict['mmr'], metric_dict['precision_recall_product'], metric_dict['clique_precision_mean'],metric_dict['clique_recall_mean'], metric_dict['clique_precision_mean_normalize'], metric_dict['clique_recall_mean_normalize'], metric_dict['clique_precision_mean_normalize_weighted'], metric_dict['clique_recall_mean_normalize_weighted'])
        #sys.stdout.flush()



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
            parameter_dict['threshold_score'] = threshold_score
            parameter_dict['inflation'] = str(inflation)
            parameter_dict['clixo_alpha'] = str(clixo_alpha)
            parameter_dict['clixo_beta'] = str(clixo_beta)
            parameter_dict['timeout'] = str(timeout)
            parameter_dict['cliquesize'] = str(cliquesize)
            parameter_dict['twostep_combination'] = twostep_combination
            parameter_dict['trim2threshold'] = trim2threshold
            parameter_dict['i'] = i
            parameter_dict['nodelete'] = nodelete 
            multiproc_input.append(parameter_dict)

        bootstrapped_cluster_predictions = p.map(cluster_helper, multiproc_input)

        #kdrew: compare full clustered set vs bootstrapped clusters
        multiproc_input = [(cluster_prediction, predicted_clusters, bootstrapped_test_networks[i]) for predicted_clusters, i in bootstrapped_cluster_predictions]
        bootstrap_cplx_cmp_metrics = p.map(comparison_helper, multiproc_input) 
        for boot_cmp in bootstrap_cplx_cmp_metrics:
            print "bootstrapped: ii %s, size %s, density %s, overlap %s, seed_method %s, fraction %s, threshold_score %s, inflation %s, clixo_alpha %s, clixo_beta %s, cliquesize %s, timeout %s, twostep_combination %s, trim2threshold: %s,  acc %s, sensitivity %s, ppv %s, mmr %s, ppi_recovered %s, precision_recall_product %s, clique_precision_mean %s, clique_recall_mean %s, clique_precision_mean_normalize %s, clique_recall_mean_normalize %s, clique_precision_mean_normalize_weighted %s, clique_recall_mean_normalize_weighted %s" % (ii, size, density, overlap, seed_method, fraction, threshold_score, inflation, clixo_alpha, clixo_beta, cliquesize, timeout, str(twostep_combination), str(trim2threshold), boot_cmp['acc'], boot_cmp['sensitivity'], boot_cmp['ppv'], boot_cmp['mmr'], boot_cmp['percent_ppi_recovered'], metric_dict['precision_recall_product'], boot_cmp['clique_precision_mean'], boot_cmp['clique_recall_mean'], metric_dict['clique_precision_mean_normalize'], metric_dict['clique_recall_mean_normalize'], metric_dict['clique_precision_mean_normalize_weighted'], metric_dict['clique_recall_mean_normalize_weighted'])


        #kdrew: keeping track of the best parameter set
        if best_eval == None or best_eval < metric_dict[args.eval_metric]: 
            best_eval = metric_dict[args.eval_metric]
            best_size = size
            best_ii = ii
            best_density = density
            best_overlap = overlap
            best_seed_method = seed_method
            best_fraction = fraction
            best_threshold_score = threshold_score
            best_inflation = inflation
            best_clixo_alpha = clixo_alpha
            best_clixo_beta = clixo_beta
            best_cliquesize = cliquesize
            best_timeout = timeout
            best_twostep_combination = twostep_combination
            best_trim2threshold = trim2threshold
            best_cluster_prediction = cluster_prediction
            print "best ii: %s size: %s density: %s overlap: %s seed_method: %s fraction: %s threshold_score: %s inflation: %s clixo_alpha: %s, clixo_beta: %s,cliquesize: %s timeout: %s twostep_combination: %s, trim2threshold: %s, numOfClusters: %s" % (best_ii, best_size, best_density, best_overlap, best_seed_method, best_fraction, best_threshold_score, best_inflation, best_clixo_alpha, best_clixo_beta, best_cliquesize, best_timeout, str(best_twostep_combination), str(trim2threshold), len(best_cluster_prediction))


        #kdrew: output best cluster prediction
        if args.output_file != None and args.output_all:
            output_filename = "%s.ii%s.%s" % (".".join(args.output_file.split('.')[:-1]), ii,  args.output_file.split('.')[-1])
            with open (output_filename, "w") as output_file:
                for cluster in cluster_prediction:
                    output_file.write(' '.join(cluster))
                    output_file.write("\n")





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
    cplx_cmp_normalize = cc.ComplexComparison(gold_std, pred_clst, normalize_by_combinations=True, pseudocount=0.00001) 

    d = dict()
    d['acc'] = cplx_cmp.acc()
    d['sensitivity'] = cplx_cmp.sensitivity()
    d['ppv'] = cplx_cmp.ppv()
    d['mmr'] = cplx_cmp.mmr()
    d['percent_ppi_recovered'] = (1.0*ppi_recovered_count) / len(test_net)
    d['precision_recall_product'] = cplx_cmp.precision_recall_product() 
    ccmm = cplx_cmp.clique_comparison_metric_mean()
    d['clique_precision_mean'] = ccmm['precision_mean']
    d['clique_recall_mean'] = ccmm['recall_mean']

    ccmm_normalize = cplx_cmp_normalize.clique_comparison_metric_mean()
    d['clique_precision_mean_normalize'] = ccmm_normalize['precision_mean']
    d['clique_recall_mean_normalize'] = ccmm_normalize['recall_mean']
    ccmm_normalize_weighted = cplx_cmp_normalize.clique_comparison_metric_mean(weighted=True)
    d['clique_precision_mean_normalize_weighted'] = ccmm_normalize_weighted['precision_mean']
    d['clique_recall_mean_normalize_weighted'] = ccmm_normalize_weighted['recall_mean']


    return d

def cluster_helper(parameter_dict):

    #print "in cluster_helper"
    #sys.stdout.flush()

    input_network_list = parameter_dict['network_list']
    ppi_scores = parameter_dict['ppi_scores']
    args = parameter_dict['args']
    size = parameter_dict['size']
    density = parameter_dict['density']
    overlap = parameter_dict['overlap']
    seed_method = parameter_dict['seed_method']
    fraction = parameter_dict['fraction']
    threshold_score = parameter_dict['threshold_score']
    i = parameter_dict['i']
    #kdrew: added some defaults, might want to think harder about the others
    clixo_alpha = parameter_dict.get('clixo_alpha',None)
    clixo_beta = parameter_dict.get('clixo_beta',None)
    inflation = parameter_dict.get('inflation',None)
    cliquesize = parameter_dict.get('cliquesize',None)
    timeout = parameter_dict.get('timeout',None)
    trim2threshold = parameter_dict.get('trim2threshold',False)
    nodelete = parameter_dict.get('nodelete',False)

    twostep_combination = parameter_dict['twostep_combination']

    #kdrew: make temp directory if does not exist
    if not os.path.exists(args.temp_dir) and args.temp_dir != None:
        os.makedirs(args.temp_dir)
    #kdrew: create temp file for bootstrapped input network, clusterone requires a file input
    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
    dirTemp = tf.mkdtemp(dir=args.temp_dir)
    try:
        #kdrew: write input_network or bootstrapped input network to temp file
        #fileTemp.write(input_network_str)
        for bstrap_ppi in input_network_list:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        if twostep_combination[0] == 'clusterone':
            #print "in clusterone"
            #sys.stdout.flush()

            #kdrew: run clustering
            proc = sp.Popen(['java', '-jar', args.clustone_jar, fileTemp.name, '-s', size, '-d', density, '--max-overlap', overlap, '--seed-method', seed_method], stdout=sp.PIPE, stderr=sp.PIPE)
            clust_out, err = proc.communicate()
            #print clust_out

            #kdrew: probably should do some error checking 
            #print err

            #kdrew: take output of cluster one (predicted clusters) and store them into list
            predicted_clusters = []
            for line in clust_out.split('\n'):
                if len(line.split() ) > 0:
                    predicted_clusters.append(line.split())

            #print "exiting clusterone"
            #sys.stdout.flush()

        elif twostep_combination[0] == 'cfinder':
            proc = sp.Popen([args.cfinder_exe, '-l', args.cfinder_license, '-i', fileTemp.name, '-o', dirTemp, '-k', cliquesize, '-t', timeout ], stdout=sp.PIPE, stderr=sp.PIPE)
            clust_out, err = proc.communicate()
            #print clust_out
            #print err
            print dirTemp

            predicted_clusters = []

            try:
                cfinder_communities_file = open("%s/k=%s/communities" % (dirTemp,cliquesize),'rb')
                for line in cfinder_communities_file.readlines():
                    if len(line.split() ) > 0 and '#' not in line:
                        #print line
                        predicted_clusters.append(line.split(':')[1].split())

            except (IOError, OSError) as e:
                print e

        elif twostep_combination[0] == 'clixo':
            proc = sp.Popen([args.clixo_bin, fileTemp.name, clixo_alpha, clixo_beta], stdout=sp.PIPE, stderr=sp.PIPE)
            clust_out, err = proc.communicate()
            #print clust_out
            #print err

            predicted_clusters = []
            for line in clust_out.split('\n'):
                if "Valid cluster:" in line:
                    gene_list = line.split()[3]
                    #print gene_list
                    predicted_clusters.append(gene_list.split(','))

    finally:
	if nodelete:
            pass
	else:
            os.remove(fileTemp.name)


    print "finished clustering phase1"
    sys.stdout.flush()

    if len(twostep_combination) >= 2:
        if twostep_combination[1] == 'mcl':
            print "MCL"
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
                                if ppi_score < threshold_score:
                                    ppi_score = 0.0

                            except KeyError:
                                ppi_score = 0.0

                            ppi_str = "%s\t%s\t%s\n" % (prot1, prot2, ppi_score)
                            fileTemp.write(ppi_str)
                        fileTemp.close()

                        proc = sp.Popen([args.mcl_bin, fileTemp.name, '--abc', '-o', outTemp.name, '-I', inflation], stdout=sp.PIPE, stderr=sp.PIPE)
                        mcl_out, err = proc.communicate()

                        #print fileTemp.name
                        #print clust
                        #print mcl_out
                        #print err

                        outfile = open(outTemp.name,"rb")
                        for line in outfile.readlines():
                        #    print line
                            mcl_clusters.append(line.split())
                        outfile.close()

                        #print "\n"


                    finally:
			if nodelete:
			    pass
			else:
			    os.remove(fileTemp.name)
                            os.remove(outTemp.name)
                        #pass
                        #print "in finally"

                predicted_clusters = mcl_clusters

        elif twostep_combination[1] == 'agglomod':
            print "AGGLOMOD"
            sys.stdout.flush()

            agglomod_clusters = []
            for clust in predicted_clusters:
                #print clust
                #sys.stdout.flush()

                graph = nx.Graph()
                #kdrew: for every pair in cluster find edge weight in input_network_list(?)
                for prot1, prot2 in it.combinations(clust,2):
                    try:
                        ppi_score = ppi_scores[frozenset([prot1,prot2])] 
                        if ppi_score < threshold_score:
                            ppi_score = 0.0

                    except KeyError:
                        ppi_score = 0.0

                    graph.add_edge(prot1,prot2,weight=ppi_score)

                #print "before NewmanGreedy"
                #sys.stdout.flush()

                newman = ag.NewmanGreedy(graph)

                #print "after NewmanGreedy"
                sys.stdout.flush()
                #print newman.quality_history
                #print newman.get_clusters()
                #sys.stdout.flush()
                agglomod_clusters += [list(x) for x in newman.get_clusters()]

            predicted_clusters  = agglomod_clusters

            #print "END OF AGGLOMOD"
            sys.stdout.flush()


    if trim2threshold:
        predicted_clusters = trim_clusters2threshold(predicted_clusters, threshold_score, ppi_scores)

    #kdrew: TODO merge complexes to elminiate redundancy
    #kdrew: use complex_merge.merge_complexes with merge_threshold 1.0

    return predicted_clusters, i


def trim_clusters2threshold(predicted_clusters, threshold_score, ppi_scores):

    trimed_clusters = []

    for clust in predicted_clusters:
        trimed_clust = []
        for prot1 in clust:
            prot1_max_score = 0.0
            for prot2 in clust:
                if prot1 != prot2:
                    try:
                        prot1_max_score = max( prot1_max_score, ppi_scores[frozenset([prot1,prot2])] )
                    except KeyError:
                        continue
            if prot1_max_score < threshold_score:
                print "removing prot1: %s max score: %s" % (prot1, prot1_max_score)
            else:
                trimed_clust.append(prot1)
        if len(trimed_clust) >  1:
            trimed_clusters.append(trimed_clust)

    return trimed_clusters




if __name__ == "__main__":
    main()


