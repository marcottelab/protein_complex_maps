
import argparse
import numpy as np
import multiprocessing as mp

import protein_complex_maps.complex_comparison as cc


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

    for i, cluster_filename in enumerate(args.cluster_filenames):

        parameter_dict = dict()
        parameter_dict['cluster_filename'] = cluster_filename
        parameter_dict['gs_boot_iteration'] = 1
        parameter_dict['gs_complexes'] = gold_standard_complexes
        parameter_dict['ex_complexes'] = excluded_complexes
        parameter_dict['id_delimin'] = args.id_delimin 
        parameter_dict['short_name'] = args.input_names[i] if args.input_names != None else None
        parameter_dict['samples'] = args.samples
        parameter_dict['exact'] = args.exact
        parameter_dict['max_clique'] = args.max_clique
        parameter_dict['pseudocount'] = args.pseudocount
        compare2goldstandard_input_list.append(parameter_dict)
    
        #compare2goldstandard(cluster_filename, gold_standard_complexes, args.id_delimin)

    #kdrew: call compare2goldstandard with pool of processors
    compare2goldstandard_results = p.map(compare2goldstandard, compare2goldstandard_input_list)

    pr_dict = dict()
    precision_dict = dict()
    recall_dict = dict()
    f1score_dict = dict()
    grand_f1score_dict = dict()
    cumulative_precision_dict = dict()
    cumulative_recall_dict = dict()
    numOfClusters_dict = dict()
    
    #kdrew: put the results into dictionaries
    for result in compare2goldstandard_results:
        precision_dict[result['ii']] = result['precision_list']
        recall_dict[result['ii']] = result['recall_list']
        f1score_dict[result['ii']] = result['f1score_list']
        grand_f1score_dict[result['ii']] = result['grand_f1score']
        cumulative_precision_dict[result['ii']] = result['cumulative_precision_list']
        cumulative_recall_dict[result['ii']] = result['cumulative_recall_list']
        numOfClusters_dict[result['ii']] = result['numOfClusters']

    #kdrew: store all dictionaries into a master dictionary
    pr_dict['precision'] = precision_dict
    pr_dict['recall'] = recall_dict
    pr_dict['f1score'] = f1score_dict
    pr_dict['grand_f1score'] = grand_f1score_dict
    pr_dict['cumulative_precision'] = cumulative_precision_dict
    pr_dict['cumulative_recall'] = cumulative_recall_dict
    pr_dict['numOfClusters'] = numOfClusters_dict

    #kdrew: output results 
    outfile = open(args.output_filename,"wb")
    outfile.write("ii$precision_list$recall_list$f1score_list$grand_f1score$cumulative_precision$cumulative_recall$numOfClusters\n" )
    for ii in pr_dict['precision'].keys():
        print ii
        precision_out = ",".join(map(str,pr_dict['precision'][ii]))
        recall_out = ",".join(map(str,pr_dict['recall'][ii]))
        f1score_out = ",".join(map(str,pr_dict['f1score'][ii]))
        grand_f1score_out = str(pr_dict['grand_f1score'][ii])
        cumulative_precision_out = ",".join(map(str,pr_dict['cumulative_precision'][ii]))
        cumulative_recall_out = ",".join(map(str,pr_dict['cumulative_recall'][ii]))
        numOfClusters_out = ",".join(map(str,pr_dict['numOfClusters'][ii]))
        outfile.write("%s$%s$%s$%s$%s$%s$%s$%s\n" % (ii,precision_out, recall_out, f1score_out, grand_f1score_out, cumulative_precision_out, cumulative_recall_out, numOfClusters_out))

    outfile.close()

#kdrew: id_delimin is for parsing a unique identifier from the filename, usually 'ii'
#def compare2goldstandard(cluster_filename, gs_complexes, id_delimin):
def compare2goldstandard(parameter_dict):
    
    cluster_filename = parameter_dict['cluster_filename']
    gs_complexes = parameter_dict['gs_complexes']
    boot_iteration = parameter_dict['gs_boot_iteration'] 
    ex_complexes = parameter_dict['ex_complexes']
    id_delimin = parameter_dict['id_delimin']
    short_name = parameter_dict['short_name']
    samples = parameter_dict['samples']
    exact = parameter_dict['exact']
    max_clique = parameter_dict['max_clique']
    pseudocount = parameter_dict['pseudocount']
    #for i, cluster_filename in enumerate(cluster_filenames):

    predicted_clusters = []
    clpred_f = open(cluster_filename,"rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())

    clpred_f.close()

    if short_name != None:
        ii = short_name
    else:
        #kdrew: strip from filename the id (ii) 
        for seg in cluster_filename.split('.'):
            if id_delimin in seg:
                ii = seg.split(id_delimin)[1]

    cplx_compare = cc.ComplexComparison(gs_complexes, predicted_clusters, exclusion_complexes=ex_complexes, samples=samples, exact=exact, max_clique=max_clique, pseudocount=pseudocount)
    
    result_dict = cplx_compare.clique_comparison_metric()
    #print result_dict
    precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
    recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]
    f1score_list = [result_dict[x]['f1score'] for x in result_dict.keys()]
    cumulative_precision_list = [result_dict[x]['cumulative_precision'] for x in result_dict.keys()]
    cumulative_recall_list = [result_dict[x]['cumulative_recall'] for x in result_dict.keys()]
    numOfClusters_list = [result_dict[x]['numOfClusters'] for x in result_dict.keys()]

    print precision_list
    print recall_list
    print f1score_list
    print cumulative_precision_list
    print cumulative_recall_list
    print numOfClusters_list
    
    #precision_dict[ii] = precision_list
    #recall_dict[ii] = recall_list
    #f1score_dict[ii] = f1score_list

    #f1score_hmean_dict[ii] = cplx_compare.clique_comparison_metric_grandf1score()
    grand_f1score = cplx_compare.clique_comparison_metric_grandf1score()
    #print f1score_hmean_dict[ii]


    #return_dict['precision'] = precision_dict
    #return_dict['recall'] = recall_dict
    #return_dict['f1score'] = f1score_dict
    #return_dict['grand_f1score'] = f1score_hmean_dict
    #
    #return return_dict

    return {'ii':ii, 'boot_iteration':boot_iteration, 'precision_list':precision_list, 'recall_list':recall_list, 'f1score_list':f1score_list, 'grand_f1score':grand_f1score, 'cumulative_precision_list':cumulative_precision_list, 'cumulative_recall_list':cumulative_recall_list, 'numOfClusters':numOfClusters_list}

#kdrew: reads in precision recall file and returns dictionary
def read_pr_file(filename, names=None):
    pr_dict1 = dict()
    pr_dict1['precision'] = dict()
    pr_dict1['recall'] = dict()
    pr_dict1['f1score'] = dict()
    pr_dict1['grand_f1score'] = dict()
    pr_dict1['cumulative_precision'] = dict()
    pr_dict1['cumulative_recall'] = dict()
    pr_dict1['numOfClusters'] = dict()

    print names

    pr_file = open(filename,"rb")
    for j, line in enumerate(pr_file.readlines()):
        if "ii$precision_list$recall_list" in line or '$' not in line:
            #kdrew: header line, do nothing
            continue
        if names != None:
            ii = names[j-1]
        else:
            ii = line.split('$')[0]
        #kdrew: fields are '$' deliminited and precision and recall lists are ',' deliminited
        precision_list = map(float,line.split('$')[1].split(','))
        recall_list = map(float,line.split('$')[2].split(','))
        f1score_list = map(float,line.split('$')[3].split(','))
        grand_f1score = float(line.split('$')[4])
        cumulative_precision_list = map(float,line.split('$')[5].split(','))
        cumulative_recall_list = map(float,line.split('$')[6].split(','))
        numOfClusters_list = map(float,line.split('$')[7].split(','))
        pr_dict1['precision'][ii] = precision_list
        pr_dict1['recall'][ii] = recall_list
        pr_dict1['f1score'][ii] = f1score_list
        pr_dict1['grand_f1score'][ii] = grand_f1score
        pr_dict1['cumulative_precision'][ii] = cumulative_precision_list
        pr_dict1['cumulative_recall'][ii] = cumulative_recall_list
        pr_dict1['numOfClusters'][ii] = numOfClusters_list

    pr_file.close()
    return pr_dict1


if __name__ == "__main__":
	main()


