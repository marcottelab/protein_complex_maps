
import argparse
import numpy as np

import protein_complex_maps.complex_comparison as cc


def main():

    parser = argparse.ArgumentParser(description="Parameter Sweep comparison to gold standard")
    parser.add_argument("--cluster_filenames", action="store", dest="cluster_filenames", nargs='+', required=True, 
                                            help="Filenames of cluster predictions")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard", required=True, 
                                            help="Filename of gold standard complexes")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Filename of where the output should go")
    args = parser.parse_args()



    #gold_standard_filename = "/home/kdrew/data/protein_complex_maps/corum_revisit/allComplexesCore_geneid_merged06_testSplit.txt"
    gold_standard_complexes = []
    gold_file = open(args.gold_standard,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())

    gold_file.close()

    pr_dict = compare2goldstandard(args.cluster_filenames, gold_standard_complexes)
    
    outfile = open(args.output_filename,"wb")
    outfile.write("ii$precision_list$recall_list\n" )
    for ii in pr_dict['precision'].keys():
        print ii
        precision_out = ",".join(map(str,pr_dict['precision'][ii]))
        recall_out = ",".join(map(str,pr_dict['recall'][ii]))
        outfile.write("%s$%s$%s\n" % (ii,precision_out, recall_out))

    outfile.close()

def compare2goldstandard(cluster_filenames, gs_complexes):
    
    #gs_plot_data = []
    return_dict = dict()
    precision_dict = dict()
    recall_dict = dict()
    
    for cluster_filename in cluster_filenames:
        predicted_clusters = []
        clpred_f = open(cluster_filename,"rb")
        for line in clpred_f.readlines():
            predicted_clusters.append(line.split())
    
        clpred_f.close()
    
        #kdrew: strip from filename the id (ii) 
        for seg in cluster_filename.split('.'):
            if 'ii' in seg:
                ii = seg.split('ii')[1]
    
        cplx_compare = cc.ComplexComparison(gs_complexes, predicted_clusters)
        
        result_dict = cplx_compare.clique_comparison_metric()
        precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
        recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]
        print precision_list
        print recall_list
        
        precision_dict[ii] = precision_list
        recall_dict[ii] = recall_list
        
        #gs_plot_data.append({'x':recall_list , 'y':precision_list, 'name':"ii: %s" % ii, 'mode' : 'lines+markers'})

    #gs_layout = Layout(title="Precision vs Recall", xaxis=XAxis(title='Recall'), yaxis=YAxis(title='Precision'))
    ##print plot_data
    #iplot({'data':gs_plot_data,'layout':gs_layout})

    return_dict['precision'] = precision_dict
    return_dict['recall'] = recall_dict

    return return_dict



if __name__ == "__main__":
	main()

