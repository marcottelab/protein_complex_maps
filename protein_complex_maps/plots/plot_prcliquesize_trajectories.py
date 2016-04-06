
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
    parser.add_argument("--cumulative_pr", action="store_true", dest="cumulative_pr", required=False, default=False,
                        help="Use cumulative precision and recall for plotting")
    parser.add_argument("--plot_numOfClusters", action="store_true", dest="plot_numOfClusters", required=False, default=False,
                        help="Plot the number of clusters that contributed to given clique size as size of plotted dot")
    args = parser.parse_args()


    pr_dict = pc.read_pr_file(args.pr_table, names=args.cluster_names)

    #colors = ['#FF7900','660000','DDB505','C2120B']

    plot_prcliquesize_trajectories(pr_dict, pr_dict['recall'].keys(), plot_filename = args.plot_filename, cumulative=args.cumulative_pr, plot_numOfClusters=args.plot_numOfClusters)



def plot_prcliquesize_trajectories(pr_dict, psweep_indices, pr_dict2=None, plot_filename=None, colors=None, cumulative=False, plot_numOfClusters=False):
    
    plot_data = []

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for ii in psweep_indices:
        numOfClusters_list = pr_dict['numOfClusters'][ii]

        if cumulative:
            precision_list = pr_dict['cumulative_precision'][ii]
            recall_list = pr_dict['cumulative_recall'][ii]
        else:
            precision_list = pr_dict['precision'][ii]
            recall_list = pr_dict['recall'][ii]

        if colors != None:
            line, = plt.plot(recall_list, precision_list, '.-', linewidth=1, label=ii, color=colors[ii])
        else:
            line, = plt.plot(recall_list, precision_list, '.-', linewidth=1, label=ii )

        for i, xy in enumerate(zip(recall_list, precision_list)):
            cliquesize = i+2
            ax.annotate('%s' % cliquesize, xy=xy, textcoords='offset points', fontsize=8) 
            if plot_numOfClusters:
                if numOfClusters_list[i] == 0:
                    markersize = 1
                else:                                                                   
                    markersize = np.log(numOfClusters_list[i])*10
                dot, = plt.plot([recall_list[i]], [precision_list[i]], '.', markersize=markersize, color=line.get_color())


    ax.set_ylim([ax.get_ylim()[0]-0.01,ax.get_ylim()[1]+0.01])
    ax.set_xlim([ax.get_xlim()[0]-0.01,ax.get_xlim()[1]+0.01])

    plt.legend(loc="best", fontsize=8)
    if cumulative:
        plt.xlabel('Cumulative Recall')
        plt.ylabel('Cumulative Precision')
    else:
        plt.xlabel('Recall')
        plt.ylabel('Precision')

    plt.savefig(plot_filename)


if __name__ == "__main__":
    main()



