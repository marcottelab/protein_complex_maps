
import argparse
import numpy as np
import multiprocessing as mp

import random
import itertools as it

import protein_complex_maps.complex_comparison as cc
import protein_complex_maps.complexes.psweep_comparison as pc



def main():

    parser = argparse.ArgumentParser(description="Parameter Sweep comparison to gold standard")
    parser.add_argument("--cluster_filenames", action="store", dest="cluster_filenames", nargs='+', required=True, 
                                            help="Filenames of cluster predictions")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard", required=True, 
                                            help="Filename of gold standard complexes")
    parser.add_argument("--excluded_complexes", action="store", dest="excluded_complexes", required=False, default=None,
                                            help="Filename of benchmark complexes to be excluded from false positive calculation")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Filename of where the output should go")
    parser.add_argument("--id_delimin", action="store", dest="id_delimin", required=False, default='ii',
                                            help="Deliminator of identifier in filename to uniquely id the predictions, default=ii (ex. predictions_ii99.txt) ")
    parser.add_argument("--input_names", action="store", nargs='+', dest="input_names", required=False, default=None, 
                                            help="Short names of input clustering, same order as cluster_filenames")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    parser.add_argument("--samples", action="store", type=int, dest="samples", required=False, default=10000,
                                    help="Number to samples for complex comparison, default=10000")
    parser.add_argument("--pseudocount", action="store", type=int, dest="pseudocount", required=False, default=1,
                                    help="Pseudocount to use during sampling (not strictly turned off for exact calculations), default=1")
    parser.add_argument("--max_clique", action="store", type=int, dest="max_clique", required=False, default=None,
                                    help="Value of largest clique size to calculate, default=None, (size of largest cluster)")
    parser.add_argument("--exact", action="store_true", dest="exact", required=False, default=False,
                                    help="Calculate the exact precision recall and f1scores rather than sampling (consider using --max_clique option as well), default=False")
    parser.add_argument("--gs_boot_iterations", action="store", type=int, dest="gs_boot_iterations", required=False, default=1,
                                    help="Number of iterations to bootstrap gold standard complexes, default=1")
    parser.add_argument("--gs_boot_fraction", action="store", type=float, dest="gs_boot_fraction", required=False, default=1.0,
                                    help="Fraction of gold standard complexes to bootstrap, default=1.0, no sampling")
    args = parser.parse_args()



    #gold_standard_filename = "/home/kdrew/data/protein_complex_maps/corum_revisit/allComplexesCore_geneid_merged06_testSplit.txt"
    gold_standard_complexes = []
    gold_file = open(args.gold_standard,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())

    gold_file.close()

    excluded_complexes = []
    if args.excluded_complexes != None:
        exclude_file = open(args.excluded_complexes,"rb")
        for line in exclude_file.readlines():
            excluded_complexes.append(line.split())

        exclude_file.close()

    p = mp.Pool(args.procs)
    #kdrew: create list of inputs to pass to parallel compare2goldstandard
    compare2goldstandard_input_list = []

    #kdrew: setup bootstrapping shuffled gold standard complexes
    sizeOfBootstrap = int(len(gold_standard_complexes)*args.gs_boot_fraction)
    #kdrew: generate k shuffled networks
    shuffled_complex_sets  = [ sorted( gold_standard_complexes, key=lambda k: random.random() ) for i in range(args.gs_boot_iterations)]
    #kdrew: split shuffled networks into partitions based on size of bootstrap
    bootstrapped_complex_sets = [ sn[:sizeOfBootstrap] for sn in shuffled_complex_sets ]

    for i, cluster_filename in enumerate(args.cluster_filenames):
        for boot_iteration, boot_gs_complexes in enumerate(bootstrapped_complex_sets):

            parameter_dict = dict()
            parameter_dict['cluster_filename'] = cluster_filename
            parameter_dict['gs_boot_iteration'] = boot_iteration
            parameter_dict['gs_complexes'] = boot_gs_complexes
            parameter_dict['ex_complexes'] = excluded_complexes
            parameter_dict['id_delimin'] = args.id_delimin 
            parameter_dict['short_name'] = args.input_names[i] if args.input_names != None else None
            parameter_dict['samples'] = args.samples
            parameter_dict['exact'] = args.exact
            parameter_dict['max_clique'] = args.max_clique
            parameter_dict['pseudocount'] = args.pseudocount
            compare2goldstandard_input_list.append(parameter_dict)
    
    #kdrew: call compare2goldstandard with pool of processors
    compare2goldstandard_results = p.map(pc.compare2goldstandard, compare2goldstandard_input_list)

    pr_dict = dict()
    precision_dict = dict()
    recall_dict = dict()
    f1score_dict = dict()
    grand_f1score_dict = dict()
    cumulative_precision_dict = dict()
    cumulative_recall_dict = dict()
    numOfClusters_dict = dict()
    precision_mean_dict = dict()
    recall_mean_dict = dict()
    precision_std_dict = dict()
    recall_std_dict = dict()
    
    #kdrew: put the results into dictionaries, where bootstrapped results are appended to a list
    for result in compare2goldstandard_results:

        try:
            precision_dict[result['ii']].append(result['precision_list'])
            recall_dict[result['ii']].append(result['recall_list'])
            numOfClusters_dict[result['ii']].append(result['numOfClusters'])
        except KeyError:
            precision_dict[result['ii']] = [result['precision_list']]
            print precision_dict[result['ii']]
            recall_dict[result['ii']] = [result['recall_list']]
            print recall_dict[result['ii']]
            numOfClusters_dict[result['ii']] = [result['numOfClusters']]


    #kdrew: for each set of predictions, combine bootstrapped samples and store mean and std
    for ii in precision_dict.keys():
        #kdrew: combines samples, leaves None type when not available
        grouped_clique_size_precisions = it.izip_longest(*precision_dict[ii])
        grouped_clique_size_recalls = it.izip_longest(*recall_dict[ii])
        precision_mean_dict[ii] = []
        precision_std_dict[ii] = []
        recall_mean_dict[ii] = []
        recall_std_dict[ii] = []
        #kdrew: for each clique size set of sampled precisions
        for clique_size_precision in grouped_clique_size_precisions:
            print clique_size_precision
            #kdrew: remove None and take mean and std
            precision_mean_dict[ii].append(np.mean([zz for zz in clique_size_precision if zz != None]))
            precision_std_dict[ii].append(np.std([zz for zz in clique_size_precision if zz != None]))

        #kdrew: for each clique size set of sampled recall
        for clique_size_recall in grouped_clique_size_recalls:
            print clique_size_recall
            #kdrew: remove None and take mean and std
            recall_mean_dict[ii].append(np.mean([zz for zz in clique_size_recall if zz != None]))
            recall_std_dict[ii].append(np.std([zz for zz in clique_size_recall if zz != None]))



    #kdrew: store all dictionaries into a master dictionary
    #pr_dict['precision'] = precision_dict
    #pr_dict['recall'] = recall_dict
    #pr_dict['f1score'] = f1score_dict
    #pr_dict['grand_f1score'] = grand_f1score_dict
    #pr_dict['cumulative_precision'] = cumulative_precision_dict
    #pr_dict['cumulative_recall'] = cumulative_recall_dict
    #pr_dict['numOfClusters'] = numOfClusters_dict

    pr_dict['precision_mean'] = precision_mean_dict
    pr_dict['precision_std'] = precision_std_dict
    pr_dict['recall_mean'] = recall_mean_dict
    pr_dict['recall_std'] = recall_std_dict

    #kdrew: output results 
    outfile = open(args.output_filename,"wb")
    outfile.write("ii$precision_mean_list$precision_std_list$recall_mean_list$recall_std_list\n" )
    for ii in pr_dict['precision_mean'].keys():
        print ii
        precision_mean_out = ",".join(map(str,pr_dict['precision_mean'][ii]))
        precision_std_out = ",".join(map(str,pr_dict['precision_std'][ii]))
        recall_mean_out = ",".join(map(str,pr_dict['recall_mean'][ii]))
        recall_std_out = ",".join(map(str,pr_dict['recall_std'][ii]))
        outfile.write("%s$%s$%s$%s$%s\n" % (ii,precision_mean_out, precision_std_out, recall_mean_out, recall_std_out,))

    outfile.close()

#kdrew: reads in precision recall file and returns dictionary
def read_pr_bootstrap_file(filename, names=None):
    pr_dict1 = dict()
    pr_dict1['precision_mean'] = dict()
    pr_dict1['precision_std'] = dict()
    pr_dict1['recall_mean'] = dict()
    pr_dict1['recall_std'] = dict()

    print names

    pr_file = open(filename,"rb")
    for j, line in enumerate(pr_file.readlines()):
        if "ii$precision_mean_list$precision_std_list" in line or '$' not in line:
            #kdrew: header line, do nothing
            continue
        if names != None:
            ii = names[j-1]
        else:
            ii = line.split('$')[0]
        #kdrew: fields are '$' deliminited and precision and recall lists are ',' deliminited
        precision_mean_list = map(float,line.split('$')[1].split(','))
        precision_std_list = map(float,line.split('$')[2].split(','))
        recall_mean_list = map(float,line.split('$')[3].split(','))
        recall_std_list = map(float,line.split('$')[4].split(','))
        pr_dict1['precision_mean'][ii] = precision_mean_list
        pr_dict1['precision_std'][ii] = precision_std_list
        pr_dict1['recall_mean'][ii] = recall_mean_list
        pr_dict1['recall_std'][ii] = recall_std_list

    pr_file.close()
    return pr_dict1


if __name__ == "__main__":
	main()


