    
import numpy as np
import pandas as pd
import pickle as p
import argparse



def main():

    parser = argparse.ArgumentParser(description="Merge complexes given a similarity threshold")
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True, 
                                            help="Filename of clusters, format one cluster per line, ids space separated")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                            help="Filename to output  merged clusters")
    parser.add_argument("--merge_threshold", action="store", type=float, dest="merge_threshold", required=True, 
                                            help="Jiccard similarity threshold on which to merge")
    parser.add_argument("--remove_largest", action="store_true", dest="remove_largest", required=False, default=False,
                                            help="Instead of merging similar clusters, remove largest")
    parser.add_argument("--complex_size_threshold", type=int, action="store", dest="complex_size", required=False, default=None,
                                    help="Size threshold for complexes, throws out complexes greater than value, default does not threshold")

    args = parser.parse_args()

    in_predicted_clusters = set()
    clpred_f = open(args.cluster_filename,"rb")
    for line in clpred_f.readlines():
        #kdrew: ignore singletons
        if len(line.split()) > 1:
            in_predicted_clusters.add(frozenset(line.split()))

    #kdrew: threshold on the number of subunits a complex can have, ie. throws out large complex, None keeps all complexes
    if args.complex_size != None:
        in_predicted_clusters = [c for c in in_predicted_clusters if len(c) <= args.complex_size]


    merged_count = 0
    final_clusters = None
    predicted_clusters = list(in_predicted_clusters)

    #kdrew: loop until there is no more clusters to merge
    while True:
        merged_clusters = set()

        print "predicted_clusters count: %s" % (len(predicted_clusters))

        was_merged = set()

        for i, cluster1 in enumerate(predicted_clusters):
            #kdrew: iterate through all clusters after cluster1
            for cluster2 in predicted_clusters[i+1:]:
                jindex = jaccard_index(cluster1,cluster2)
                if jindex >= args.merge_threshold:
                    cluster1_set = frozenset(cluster1)
                    cluster2_set = frozenset(cluster2)
                    #print "merging %s : %s" % (cluster1, cluster2)

                    if args.remove_largest:
                        #kdrew: add the smaller cluster to the final set
                        if len(cluster1_set) < len(cluster2_set):
                            merged_clusters.add(cluster1_set)
                        else:
                            merged_clusters.add(cluster2_set)

                    else:
                        cluster_union = frozenset(cluster1_set.union(cluster2_set))
                        merged_clusters.add(cluster_union)

                    merged_count += 1

                    #kdrew: remember the clusters that were merged
                    was_merged.add(cluster1_set)
                    was_merged.add(cluster2_set)

                    break

        print "was_merged count: %s" % (len(was_merged))
        for cluster in predicted_clusters:
            #kdrew: if wasn't merged add to set
            if frozenset(cluster) not in was_merged:
                merged_clusters.add(frozenset(cluster))


        #kdrew: test for convergence and reset counter
        print "merged_count: %s" % (merged_count)
        if merged_count == 0:
            final_clusters = list(merged_clusters)
            break
        else:
            merged_count = 0
            predicted_clusters = list(merged_clusters)

    outfile = open(args.output_filename,"wb")
    for cluster in final_clusters:
        outfile.write(" ".join(cluster) + '\n')

    outfile.close()
                

def jaccard_index(x, y):
    sx = set(x)
    sy = set(y)
    return 1.0 * len(sx.intersection(sy)) / len(sx.union(sy))



if __name__ == "__main__":
        main()

