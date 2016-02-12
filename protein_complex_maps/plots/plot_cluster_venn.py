
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2

import itertools as it

import argparse

def main():

    parser = argparse.ArgumentParser(description="Plot Venn diagram of input clusterings")

    parser.add_argument("--input_cluster_files", action="store", nargs='+', dest="cluster_files", required=True, 
                        help="Names of clustering files, at most 3")
    parser.add_argument("--input_cluster_names", action="store", nargs='+', dest="cluster_names", required=False, default=None, 
                        help="Names of clustering, same order as cluster_files, to be used for plotting legends")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                        help="Filename of output plot")
    parser.add_argument("--clique_size", action="store", dest="clique_size", type=int, required=False, default=2,
                        help="Clique size to use for overlap, default=2")
    args = parser.parse_args()


    fnames = args.cluster_files
    cnames = args.cluster_names

    clusterings_dict = dict()
    for filename in args.cluster_files:
        clusters = []
        f = open(filename,"rb")
        for line in f.readlines():
            clusters.append(set(line.split()))

        clusterings_dict[filename] = clusters

    clique_dict = dict()
    for k in clusterings_dict.keys():
        clique_set = set()
        for clust in clusterings_dict[k]:
            clust_clique_set = [x for x in it.combinations(clust,args.clique_size)]
            clique_set = clique_set.union(clust_clique_set)
        clique_dict[k] = clique_set

    if len(args.cluster_files) == 2:
        c0 = clique_dict[fnames[0]]
        c1 = clique_dict[fnames[1]]
        both_shared = c0.intersection(c1)
        both_count = len(both_shared)
        c0_count = len(c0) - both_count
        c1_count = len(c1) - both_count
        venn2(subsets = (c0_count, c1_count, both_count), set_labels = (fnames[0],fnames[1]))
    
    if len(args.cluster_files) == 3:
        c0 = clique_dict[fnames[0]]
        c1 = clique_dict[fnames[1]]
        c2 = clique_dict[fnames[2]]

        all_shared = c0.intersection(c1.intersection(c2)) 

        c0_c1_shared = c0.intersection(c1)
        c0_c2_shared = c0.intersection(c2)
        c1_c2_shared = c1.intersection(c2)

        print c0_c2_shared - all_shared
    
        all_count = len(all_shared)
        c0_c1_count = len(c0_c1_shared) - all_count
        c0_c2_count = len(c0_c2_shared) - all_count
        c1_c2_count = len(c1_c2_shared) - all_count
        c0_count = len(c0) - c0_c1_count - c0_c2_count - all_count
        c1_count = len(c1) - c0_c1_count - c1_c2_count - all_count
        c2_count = len(c2) - c0_c2_count - c1_c2_count - all_count
        
        venn3(subsets = (c0_count, c1_count, c0_c1_count, c2_count, c0_c2_count, c1_c2_count, all_count), set_labels = (cnames[0],cnames[1],cnames[2]))

    #all_only = 291
    #hein_bioplex_only = 927 - all_only
    #hein_fraction_only = 3078 - all_only
    #fraction_bioplex_only = 2954 - all_only
    #hein_only = 25543 - hein_bioplex_only - hein_fraction_only - all_only
    #bioplex_only = 52366 - hein_bioplex_only - fraction_bioplex_only - all_only
    #fraction_only = 999424 - hein_fraction_only - fraction_bioplex_only - all_only
    #venn3(subsets = (fraction_only, hein_only, hein_fraction_only, bioplex_only, fraction_bioplex_only, hein_bioplex_only, all_only), set_labels = ('Fractionation', 'Hein', 'Bioplex'))


    plt.savefig(args.plot_filename)


if __name__ == "__main__":
    main()



