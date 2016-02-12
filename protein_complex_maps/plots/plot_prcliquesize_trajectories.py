
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import gridspec

#plt.style.use('ggplot')

import protein_complex_maps.complex_comparison as cc
import protein_complex_maps.complexes.psweep_comparison as pc

import numpy as np
import statsmodels.api as sm

from scipy.stats import hmean
from scipy.stats import spearmanr

import operator

import argparse

def main():

    parser = argparse.ArgumentParser(description="Plot F1-score vs clique size")

    parser.add_argument("--input_pr_table", action="store", dest="pr_table", required=True, 
                        help="Names of precision recall table file")
    parser.add_argument("--input_cluster_names", action="store", nargs='+', dest="cluster_names", required=False, default=None, 
                        help="Names of clustering, same order as pr_table, to be used for plotting legends")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                        help="Filename of output plot")
    args = parser.parse_args()


    pr_dict = pc.read_pr_file(args.pr_table, names=args.cluster_names)

    plot_prcliquesize_trajectories(pr_dict, pr_dict['recall'].keys(), plot_filename = args.plot_filename)



def plot_prcliquesize_trajectories(pr_dict, psweep_indices, pr_dict2=None, plot_filename=None):
    
    plot_data = []

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for ii in psweep_indices:
        precision_list = pr_dict['precision'][ii]
        recall_list = pr_dict['recall'][ii]

        line, = plt.plot(recall_list, precision_list, '.-', linewidth=1, label=ii)

        for i, xy in enumerate(zip(recall_list, precision_list)):
            cliquesize = i+2
            ax.annotate('%s' % cliquesize, xy=xy, textcoords='offset points', fontsize=6) 


    plt.legend(loc="lower right", fontsize=8)
    plt.xlabel('Recall')
    plt.ylabel('Precision')

    plt.savefig(plot_filename)


if __name__ == "__main__":
    main()



