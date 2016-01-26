
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
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Filename of where the output should go")
    parser.add_argument("--id_delimin", action="store", dest="id_delimin", required=False, default='ii',
                                            help="Deliminator of identifier in filename to uniquely id the predictions, default=ii (ex. predictions_ii99.txt) ")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    parser.add_argument("--samples", action="store", type=int, dest="samples", required=False, default=10000,
                                    help="Number to samples for complex comparison, defaul=10000")
    args = parser.parse_args()



    #gold_standard_filename = "/home/kdrew/data/protein_complex_maps/corum_revisit/allComplexesCore_geneid_merged06_testSplit.txt"
    gold_standard_complexes = []
    gold_file = open(args.gold_standard,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())

    gold_file.close()

    p = mp.Pool(args.procs)
    #kdrew: create list of inputs to pass to parallel compare2goldstandard
    compare2goldstandard_input_list = []
    for i, cluster_filename in enumerate(args.cluster_filenames):
        parameter_dict = dict()
        parameter_dict['cluster_filename'] = cluster_filename
        parameter_dict['gs_complexes'] = gold_standard_complexes
        parameter_dict['id_delimin'] = args.id_delimin
        parameter_dict['samples'] = args.samples
        compare2goldstandard_input_list.append(parameter_dict)
    
        #compare2goldstandard(cluster_filename, gold_standard_complexes, args.id_delimin)

    #kdrew: call compare2goldstandard with pool of processors
    compare2goldstandard_results = p.map(compare2goldstandard, compare2goldstandard_input_list)

    pr_dict = dict()
    precision_dict = dict()
    recall_dict = dict()
    f1score_dict = dict()
    grand_f1score_dict = dict()
    
    #kdrew: put the results into dictionaries
    for result in compare2goldstandard_results:
        precision_dict[result['ii']] = result['precision_list']
        recall_dict[result['ii']] = result['recall_list']
        f1score_dict[result['ii']] = result['f1score_list']
        grand_f1score_dict[result['ii']] = result['grand_f1score']

    #kdrew: store all dictionaries into a master dictionary
    pr_dict['precision'] = precision_dict
    pr_dict['recall'] = recall_dict
    pr_dict['f1score'] = f1score_dict
    pr_dict['grand_f1score'] = grand_f1score_dict

    #kdrew: output results 
    outfile = open(args.output_filename,"wb")
    outfile.write("ii$precision_list$recall_list$f1score_list$grand_f1score\n" )
    for ii in pr_dict['precision'].keys():
        print ii
        precision_out = ",".join(map(str,pr_dict['precision'][ii]))
        recall_out = ",".join(map(str,pr_dict['recall'][ii]))
        f1score_out = ",".join(map(str,pr_dict['f1score'][ii]))
        grand_f1score_out = str(pr_dict['grand_f1score'][ii])
        outfile.write("%s$%s$%s$%s$%s\n" % (ii,precision_out, recall_out, f1score_out, grand_f1score_out))

    outfile.close()

#kdrew: id_delimin is for parsing a unique identifier from the filename, usually 'ii'
#def compare2goldstandard(cluster_filename, gs_complexes, id_delimin):
def compare2goldstandard(parameter_dict):
    
    cluster_filename = parameter_dict['cluster_filename']
    gs_complexes = parameter_dict['gs_complexes']
    id_delimin = parameter_dict['id_delimin']
    samples = parameter_dict['samples']
    #for i, cluster_filename in enumerate(cluster_filenames):

    predicted_clusters = []
    clpred_f = open(cluster_filename,"rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())

    clpred_f.close()

    #kdrew: strip from filename the id (ii) 
    for seg in cluster_filename.split('.'):
        if id_delimin in seg:
            ii = seg.split(id_delimin)[1]

    cplx_compare = cc.ComplexComparison(gs_complexes, predicted_clusters, samples=samples)
    
    result_dict = cplx_compare.clique_comparison_metric()
    #print result_dict
    precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
    recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]
    f1score_list = [result_dict[x]['f1score'] for x in result_dict.keys()]

    print precision_list
    print recall_list
    print f1score_list
    
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

    return {'ii':ii, 'precision_list':precision_list, 'recall_list':recall_list, 'f1score_list':f1score_list, 'grand_f1score':grand_f1score}

#kdrew: reads in precision recall file and returns dictionary
def read_pr_file(filename):
    pr_dict1 = dict()
    pr_dict1['precision'] = dict()
    pr_dict1['recall'] = dict()
    pr_dict1['f1score'] = dict()
    pr_dict1['grand_f1score'] = dict()

    pr_file = open(filename,"rb")
    for line in pr_file.readlines():
        if "ii$precision_list$recall_list" in line:
            #kdrew: header line, do nothing
            continue
        ii = line.split('$')[0]
        #kdrew: fields are '$' deliminited and precision and recall lists are ',' deliminited
        precision_list = map(float,line.split('$')[1].split(','))
        recall_list = map(float,line.split('$')[2].split(','))
        f1score_list = map(float,line.split('$')[3].split(','))
        grand_f1score = float(line.split('$')[4])
        pr_dict1['precision'][ii] = precision_list
        pr_dict1['recall'][ii] = recall_list
        pr_dict1['f1score'][ii] = f1score_list
        pr_dict1['grand_f1score'][ii] = grand_f1score

    pr_file.close()
    return pr_dict1


if __name__ == "__main__":
	main()


