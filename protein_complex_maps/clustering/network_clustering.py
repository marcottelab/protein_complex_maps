
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

import networkx as nx
import agglomcluster.agglomod as ag

import protein_complex_maps.evaluation.complex_comparison as cc

#kdrew: original code from clustering_parameter_optimization.py
#kdrew: this script just does clustering on the input network for the input parameters, no evaluation

def main():

    parser = argparse.ArgumentParser(description="Finds optimal parameters for clustering of ppi network")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True, 
                                    help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output clusters for best set of parameters")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", required=False, default=None,
                                    help="Sets random seed (int), default=None")

    parser.add_argument("--ppi_threshold_score", action="store", dest="ppi_threshold_score", nargs='+', required=False, 
                                    default=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],
                                    help="Use ppis with score higher or equal to threshold score, default = 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1")

    parser.add_argument("--clusterone_jar", action="store", dest="clustone_jar", required=False, 
                                    default="/home/kdrew/programs/clusterone/cluster_one-1.0.jar",
                                    help="Location of cluster_one jar file, default= /home/kdrew/programs/clusterone/cluster_one-1.0.jar")
    parser.add_argument("--clusterone_size", action="store", dest="clusterone_size", nargs='+', required=False, 
                                    default=[2,],
                                    help="ClusterOne Size parameter sweep, default = 2")
    parser.add_argument("--clusterone_density", action="store", dest="clusterone_density", nargs='+', required=False, 
                                    default=[.25], 
                                    help="ClusterOne Density parameter sweep, default = .25")
    parser.add_argument("--clusterone_max_overlap", action="store", dest="clusterone_max_overlap", nargs='+', required=False, 
                                    default=[0.75], 
                                    help="ClusterOne max-overlap parameter sweep, default = 0.75")
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

    parser.add_argument("--mcl_bin", action="store", dest="mcl_bin", required=False, 
                                    default='/usr/bin/mcl',
                                    help="Location of mcl binary, default = '/usr/bin/mcl' ")
    parser.add_argument("--mcl_inflation", action="store", dest="mcl_inflation", nargs='+', required=False, 
                                    default=[3],
                                    help="MCL Inflation (-I) parameter, default = 3")

    parser.add_argument("--cfinder_bin", action="store", dest="cfinder_bin", required=False, 
                                    default="/home/kdrew/programs/CFinder-2.0.6--1448/CFinder_commandline64",
                                    help="Location of CFinder executable, default= /home/kdrew/programs/CFinder-2.0.6--1448/CFinder_commandline64")
    parser.add_argument("--cfinder_license", action="store", dest="cfinder_license", required=False, 
                                    default="/home/kdrew/programs/CFinder-2.0.6--1448/licence.txt",
                                    help="Location of CFinder license, default= /home/kdrew/programs/CFinder-2.0.6--1448/licence.txt")
    parser.add_argument("--cfinder_cliquesize", action="store", dest="cfinder_cliquesize", nargs='+', required=False, 
                                    default=[3],
                                    help="Cfinder clique size (-k) parameter, default = 3 (CFinder's default setting)")
    parser.add_argument("--cfinder_timeout", action="store", dest="cfinder_timeout", nargs='+', required=False, 
                                    default=[10],
                                    help="Cfinder timeout (-t) parameter, default = 10 (CFinder's default setting)")

    parser.add_argument("--trim2threshold", action="store_true", dest="trim2threshold", required=False, default=True,
                                    help="Trim final clusters of subunits that do not have an edge that passes the threshold_score, default=True")
    parser.add_argument("--twostep_combination", action="store", dest="twostep_combination", nargs='+', required=False, 
                                    default=['clusterone','mcl'],
                                    help="Combination of two step clustering, default = [clusterone,mcl], options=[clusterone,mcl,cfinder,agglomod,clixo]")
    parser.add_argument("--unweighted_one_step", action="store_true", dest="unweighted_one_step", required=False, default=False,
                                    help="Do not use weights for the first step of clustering, default = False")

    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    parser.add_argument("--temp_dir", action="store", dest="temp_dir", required=False, default=None,
                                    help="Where to store temporary files generated during processing, default=None (defaults to OS level tmp), in memory suggestion = /dev/shm/")
    parser.add_argument("--nodelete", action="store_true", dest="nodelete", required=False, default=False,
                                    help="When set, does not delete temporary files")
    parser.add_argument("--starting_id_num", action="store", type=int, dest="starting_id_num", required=False, default=0,
                                    help="Start the clustering identifier (ii) at a specific integer, default=0)")


    args = parser.parse_args()

    if 'clusterone' in args.twostep_combination:
        if not os.path.isfile(args.clustone_jar):
            raise IOError("File not found: %s" % args.clustone_jar)

    if 'mcl' in args.twostep_combination:
        if not os.path.isfile(args.mcl_bin):
            raise IOError("File not found: %s" % args.mcl_bin)

    if 'cfinder' in args.twostep_combination:
        if not os.path.isfile(args.cfinder_bin):
            raise IOError("File not found: %s" % args.cfinder_bin)
        if not os.path.isfile(args.cfinder_license):
            raise IOError("File not found: %s" % args.cfinder_licence)

    if 'clixo' in args.twostep_combination:
        if not os.path.isfile(args.clixo_bin):
            raise IOError("File not found: %s" % args.clixo_bin)
        if args.clixo_alpha == None:
            print "WARNING: --clixo_alpha not set"
        if args.clixo_beta == None:
            print "WARNING: --clixo_beta not set"

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

    p = mp.Pool(args.procs)
    network_input_list = []
    #kdrew: enumerate through all combinations of inputs
    for ii, parameters  in enumerate(it.product(args.clusterone_size, args.clusterone_density, args.clusterone_max_overlap, args.clusterone_seed_method, args.mcl_inflation, args.clixo_alpha, args.clixo_beta, args.cfinder_cliquesize, args.cfinder_timeout, args.ppi_threshold_score)):

        #kdrew: start ii at a different number other than default 0
        ii = ii+args.starting_id_num

        #/print ii, parameters
        sys.stdout.flush()
        #kdrew: unpack parameters
        size, density, overlap, seed_method, inflation, clixo_alpha, clixo_beta, cliquesize, timeout, threshold_score = parameters
        threshold_score = float(threshold_score)

        #print "using threshold_score: %s" % (len(input_network_list))
        network_list = []
        for x in input_network_list:                   
            if float(x.split()[2]) >= float(threshold_score):
                network_list.append(x)
            else:
                #kdrew: assuming sorted
                break

        fraction = 1.0*len(network_list) / len(input_network_list)

        #kdrew: setting up parameters to be sent to different threads, many of the parameters end up on commandlines so they are converted to strings here
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
        parameter_dict['unweighted_one_step'] = args.unweighted_one_step
        parameter_dict['ii'] = ii
        parameter_dict['nodelete'] = args.nodelete 

        #kdrew: append to list of parameters for input into multiprocessor map
        network_input_list.append(parameter_dict) 

    #kdrew: write out parameter file
    cols = ['ii', 'size', 'density', 'overlap', 'seed_method', 'fraction', 'threshold_score',  'inflation', 'clixo_alpha', 'clixo_beta', 'cliquesize', 'timeout', 'twostep_combination', 'trim2threshold','unweighted_one_step']
    ni_df = pd.DataFrame(network_input_list)[cols]
    parameter_filename = "%s.params" % (".".join(args.output_file.split('.')[:-1]))
    ni_df.to_csv(parameter_filename, index=False)

    #kdrew: call clustering on parameter combinations
    p.map(cluster_helper, network_input_list)

    p.close()
    p.join()

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
    ii = parameter_dict['ii']
    #kdrew: added some defaults, might want to think harder about the others
    clixo_alpha = parameter_dict.get('clixo_alpha',None)
    clixo_beta = parameter_dict.get('clixo_beta',None)
    inflation = parameter_dict.get('inflation',None)
    cliquesize = parameter_dict.get('cliquesize',None)
    timeout = parameter_dict.get('timeout',None)
    trim2threshold = parameter_dict.get('trim2threshold',False)
    unweighted_one_step = parameter_dict.get('unweighted_one_step',False)

    nodelete = parameter_dict.get('nodelete',False)

    twostep_combination = parameter_dict['twostep_combination']

    #kdrew: make temp directory if does not exist
    if not os.path.exists(args.temp_dir) and args.temp_dir != None:
        os.makedirs(args.temp_dir)
    #kdrew: create temp file for input network, clusterone requires a file input
    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
    dirTemp = tf.mkdtemp(dir=args.temp_dir)
    try:
        #kdrew: write input_network to temp file
        #fileTemp.write(input_network_str)
        for ppi in input_network_list:
            if unweighted_one_step:
                fileTemp.write("%s\n" % ' '.join(ppi.split()[:2]))
            else:
                fileTemp.write(ppi)
        fileTemp.close()

        if twostep_combination[0] == 'clusterone':
            #print "in clusterone"
            #sys.stdout.flush()

            #kdrew: run clustering
            proc = sp.Popen(['java', '-jar', args.clustone_jar, fileTemp.name, '-s', size, '-d', density, '--max-overlap', overlap, '--seed-method', seed_method], stdout=sp.PIPE, stderr=sp.PIPE)
            clust_out, err = proc.communicate()
            print clust_out

            #kdrew: probably should do some error checking 
            print err

            #kdrew: take output of cluster one (predicted clusters) and store them into list
            predicted_clusters = []
            for line in clust_out.split('\n'):
                if len(line.split() ) > 0:
                    predicted_clusters.append(line.split())

            #print "exiting clusterone"
            #sys.stdout.flush()

        elif twostep_combination[0] == 'cfinder':
            proc = sp.Popen([args.cfinder_bin, '-l', args.cfinder_license, '-i', fileTemp.name, '-o', dirTemp, '-k', cliquesize, '-t', timeout ], stdout=sp.PIPE, stderr=sp.PIPE)
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
            mcl_clusters = []
            #kdrew: for each predicted cluster, recluster using MCL
            for clust in predicted_clusters:
                #kdrew: for every pair in cluster find edge weight in input_network_list(?)
                fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
                outTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)

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

    #kdrew: output clusters
    if args.output_file != None:
        output_filename = "%s.ii%s.%s" % (".".join(args.output_file.split('.')[:-1]), ii,  args.output_file.split('.')[-1])
        with open (output_filename, "w") as output_file:
            for cluster in predicted_clusters:
                output_file.write(' '.join(cluster))
                output_file.write("\n")

    #return predicted_clusters, ii


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


