
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import gridspec

#plt.style.use('ggplot')

import protein_complex_maps.complexes.psweep_comparison as pc
import protein_complex_maps.plots.plot_prcliquesize_trajectories as ppt

import argparse

def main():

    parser = argparse.ArgumentParser(description="Plot F1-score vs clique size with specified color scheme")

    parser.add_argument("--input_pr_table", action="store", dest="pr_table", required=True, 
                        help="Names of precision recall table file")
    parser.add_argument("--input_cluster_names", action="store", nargs='+', dest="cluster_names", required=False, default=None, 
                        help="Names of clustering, same order as pr_table, to be used for plotting legends")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                        help="Filename of output plot")
    parser.add_argument("--colors", action="store", nargs='+', dest="colors", required=False, default=['#C2120B','#660000','#DDB505','#FF7900'], 
                        help="Colors of plot lines, same order as pr_table, to be used for plotting legends")
    args = parser.parse_args()


    pr_dict = pc.read_pr_file(args.pr_table, names=args.cluster_names)

    color_dict = dict()
    for i, name in enumerate(args.cluster_names):
        color_dict[name] = args.colors[ i % len(args.colors) ]

    ppt.plot_prcliquesize_trajectories(pr_dict, pr_dict['recall'].keys(), plot_filename = args.plot_filename, colors=color_dict)


if __name__ == "__main__":
    main()



