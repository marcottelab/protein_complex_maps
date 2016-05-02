
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import gridspec

import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

#plt.style.use('ggplot')

import protein_complex_maps.complex_comparison as cc
import protein_complex_maps.complexes.psweep_comparison as pc
import protein_complex_maps.complexes.psweep_comparison_bootstrap as pcb

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
    parser.add_argument("--bootstrap", action="store_true", dest="bootstrap", required=False, default=False,
                        help="Plot bootstrapped precision recall mean and std, requires input_pr_table run from psweep_comparison_bootstrap.py")
    parser.add_argument("--errorbars", action="store_true", dest="errorbars", required=False, default=False,
                        help="Plot errorbars for bootstrapped precision recall mean and std, requires input_pr_table run from psweep_comparison_bootstrap.py")
    parser.add_argument("--colors", action="store", nargs='+', dest="colors", required=False, default=None,
                        help="Colors of plot lines, same order as pr_table, to be used for plotting legends, example: #C2120B #660000 #DDB505 #FF7900")
    args = parser.parse_args()


    color_dict = None
    if args.colors != None:
        color_dict = dict()
        for i, name in enumerate(args.cluster_names):
            color_dict[name] = args.colors[ i % len(args.colors) ]
    print color_dict

    if args.bootstrap:
        pr_dict = pcb.read_pr_bootstrap_file(args.pr_table, names=args.cluster_names)
        plot_prcliquesize_trajectories(pr_dict, pr_dict['recall_mean'].keys(), plot_filename = args.plot_filename, bootstrap=args.bootstrap, colors=color_dict, errorbars=args.errorbars)
    else:
        pr_dict = pc.read_pr_file(args.pr_table, names=args.cluster_names)
        plot_prcliquesize_trajectories(pr_dict, pr_dict['recall'].keys(), plot_filename = args.plot_filename, cumulative=args.cumulative_pr, plot_numOfClusters=args.plot_numOfClusters, colors=color_dict)

    #colors = ['#FF7900','660000','DDB505','C2120B']


def plot_prcliquesize_trajectories(pr_dict, psweep_indices, pr_dict2=None, plot_filename=None, colors=None, cumulative=False, plot_numOfClusters=False, bootstrap=False, errorbars=False):
    
    plot_data = []

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for ii in psweep_indices:
        if plot_numOfClusters:
            numOfClusters_list = pr_dict['numOfClusters'][ii]

        if cumulative:
            precision_list = pr_dict['cumulative_precision'][ii]
            recall_list = pr_dict['cumulative_recall'][ii]
        elif bootstrap:
            precision_list = pr_dict['precision_mean'][ii]
            recall_list = pr_dict['recall_mean'][ii]
            precision_std_list = pr_dict['precision_std'][ii]
            recall_std_list = pr_dict['recall_std'][ii]
            xerr=None
            yerr=None
            if errorbars:
                xerr = recall_std_list
                yerr = precision_std_list
            if colors != None:
                ax.errorbar(recall_list, precision_list, xerr=xerr, yerr=yerr, fmt='-o', label=ii, markersize=4, color=colors[ii])
            else:
                ax.errorbar(recall_list, precision_list, xerr=xerr, yerr=yerr, fmt='-o', markersize=4, label=ii)

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
                #if numOfClusters_list[i] == 0:
                #    markersize = 1
                #else:                                                                   
                #    markersize = np.log(numOfClusters_list[i])*10
                #dot, = plt.plot([recall_list[i]], [precision_list[i]], '.', markersize=markersize, color=line.get_color())

                radius = np.sqrt(1.0*numOfClusters_list[i] / np.pi)

                #kdrew: from http://stackoverflow.com/questions/9230389/why-is-matplotlib-plotting-my-circles-as-ovals
                # calculate asymmetry of x and y axes:
                x0, x1 = ax.get_xlim()
                y0, y1 = ax.get_ylim()
                dx = x1 - x0
                dy = y1 - y0
                #width = maxd / dx * radius/250
                maxd = max(dx,dy)
                width =  dx / maxd * radius/50
                height = dy / maxd * radius/50
                
                #circle = plt.Circle((recall_list[i], precision_list[i]), radius/250, color=line.get_color(), fill=False)
                #fig = plt.gcf()
                #fig.gca().add_artist(circle)

                #circle = mpatches.Circle((recall_list[i], precision_list[i]), radius/250, ec="none", transform=None)
                circle = mpatches.Ellipse((recall_list[i], precision_list[i]), width, height, ec="none", transform=None)
                collection = PatchCollection([circle], color=line.get_color(), alpha=0.3)
                ax.add_collection(collection)




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



