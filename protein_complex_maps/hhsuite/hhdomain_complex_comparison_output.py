
import argparse
import numpy as np
import operator
import pickle as p

import protein_complex_maps.hhsuite.hhdomain as hhd

INF_ZSCORE = 100.0

def main():

    parser = argparse.ArgumentParser(description="Filter and output complex comparison tables")
    parser.add_argument("--results_pickle_filename", action="store", dest="results_pickle_filename", required=True, 
                                            help="Filename to read results in pickle file")
    parser.add_argument("--complex_filename", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                            help="Filename to store scoring output between clusters")
    parser.add_argument("--protein_edge_output_filename", action="store", dest="protein_edge_output_filename", required=True, 
                                            help="Filename to store output protein edges")
    parser.add_argument("--clustering_filename", action="store", dest="clustering_filename", required=False, default=None, 
                                            help="Filename of clusters of complexes, limits edges to be output between complexes in the same cluster, (not implemented yet)")
    parser.add_argument("--overlap_score_thershold", action="store", type=int, dest="overlap_score_thershold", required=False, default=3,
                                            help="Output results have an overlap score (number of common domains) atleast threshold, default=3")
    parser.add_argument("--protein_count_threshold", action="store", type=int, dest="protein_count_threshold", required=False, default=1,
                                            help="Output results have an protein count (number of proteins with common domains) atleast threshold, default=1")
    parser.add_argument("--jaccard_threshold", action="store", type=float, dest="jaccard_threshold", required=False, default=1.0,
                                            help="Output results have jaccard index less than threshold between complexes, default=1.0 (no thresholding)")
    parser.add_argument("--sep", action="store", dest="sep", type=str, required=False, default=',',
                                            help="Separator for complex comparison results, default=','")
    parser.add_argument("--output_fields", action="store", dest="output_fields", nargs='+', required=False, default=['overlap_score','zscore','jaccard','protein_count'],
                                            help="Specify the fields to output (overlap_score, zscore, jaccard, protein_count) default=all")
    args = parser.parse_args()

    complexes = []
    complex_file = open(args.complex_filename,"rb")
    for line in complex_file.readlines():
        complexes.append(line.split())
    complex_file.close()

    comparison_results = p.load(open(args.results_pickle_filename,"rb"))

    fout = open(args.output_filename,"wb")
    prot_out = open(args.protein_edge_output_filename,"wb")
    #fout.write("i, j, overlap_score, zscore, jaccard, protein_count\n") 
    header_list = ['i','j'] + args.output_fields
    fout.write("%s\n" % args.sep.join(header_list)) 
    for i,j in comparison_results.keys():
        jaccard_coefficient = jaccard_index(complexes[i],complexes[j]) 
        if comparison_results[(i,j)]['overlap_score'] >= args.overlap_score_thershold and comparison_results[(i,j)]['protein_count'] >= args.protein_count_threshold and jaccard_coefficient <= args.jaccard_threshold:

            #fout.write("%s,%s,%s,%s,%s,%s\n" % (i,j,comparison_results[(i,j)]['overlap_score'],
            #                                                comparison_results[(i,j)]['zscore'], 
            #                                                jaccard_coefficient,
            #                                                comparison_results[(i,j)]['protein_count'],))

            out_list = [i,j]
            for field in args.output_fields:
                if 'overlap_score' == field: 
                    out_list.append(comparison_results[(i,j)]['overlap_score'])
                if 'zscore' == field:
                    out_list.append(comparison_results[(i,j)]['zscore'])
                if 'jaccard' == field:
                    out_list.append(jaccard_coefficient)
                if 'protein_count' == field:
                    out_list.append(comparison_results[(i,j)]['protein_count'])

            out_str = args.sep.join(map(str,out_list))
            fout.write( "%s\n" % out_str )

            for domain_pair in comparison_results[(i,j)]['overlap_domains']:
                #kdrew: output is in cytoscape edge format
                prot_out.write("%s_%s\t(%s_%s)\t%s_%s\t%s\n" % (i, domain_pair[0][0], domain_pair[0][1], domain_pair[0][2], j, domain_pair[1][0], domain_pair[0][0] == domain_pair[1][0]))

    fout.close()
    prot_out.close()



def jaccard_index(x, y):
    sx = set(x)
    sy = set(y)
    return 1.0 * len(sx.intersection(sy)) / len(sx.union(sy))


if __name__ == "__main__":
	main()


