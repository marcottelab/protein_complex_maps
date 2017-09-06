
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
    #parser.add_argument("--clique_size", action="store", dest="clique_size", type=int, required=False, default=2,
    #                    help="Clique size to use for overlap, default=2")
    args = parser.parse_args()


    pr_dict = pc.read_pr_file(args.pr_table, names=args.cluster_names)

    plotF1score_vs_cliqueSize(pr_dict, pr_dict['recall'].keys(), plot_filename = args.plot_filename)



def plotF1score_vs_cliqueSize(pr_dict, psweep_indices, pr_dict2=None, plot_filename=None):
    
    plot_data = []

    #gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    #fig = plt.figure()
    #ax = fig.add_subplot(2,1,1, adjustable='box', aspect=0.8)
    #ax = fig.add_subplot(1,2,1, adjustable='box', aspect=0.8)
    #ax2 = fig.add_subplot(2,1,2, adjustable='box', aspect=0.2)
    #ax2 = fig.add_subplot(1,2,2, adjustable='box')
    #ax2= plt.subplots(gs[1], sharey=ax)
    #f, (ax, ax2) = plt.subplots(1, 2, sharey=True)


    #ax1 = plt.subplot2grid((3, 3), (0, 0))
    #ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=2)
    #ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    #ax4 = plt.subplot2grid((3, 3), (1, 2), rowspan=2)

    ax = plt.subplot2grid((1, 10), (0, 0), colspan=9)
    ax2 = plt.subplot2grid((1, 10), (0, 9))



    for ii in psweep_indices:
        f1score_list = pr_dict['f1score'][ii]
        clique_list = range(2,len(f1score_list)+2)
        #plot_data.append({'x':clique_list , 'y':f1score_list, 'name':"%s" % ii.split('/')[-1], 'mode' : 'lines+markers', 'line' : line})

        line, = ax.plot(clique_list, f1score_list, '.-', linewidth=1, label=ii)
        line, = ax2.plot(clique_list, f1score_list, '.-', linewidth=1, label=ii)

    #layout = Layout(title="F1-Score vs Clique Size", xaxis=XAxis(title='Clique Size'), yaxis=YAxis(title='F1-Score'), legend=dict(x=0.4,y=1))
    #print plot_data
    #iplot({'data':plot_data,'layout':layout})
    #if plot_filename != None:
    #    py.image.save_as({'data': plot_data, 'layout':layout}, plot_filename)
 

    #kdrew: hardcoded for published complex + ii149 vs corum testsplit1
    ax.set_ylim([-0.01,0.4])
    ax.set_xlim([-0.01,22])
    ax2.set_ylim([-0.01,0.4])
    ax2.set_xlim([38,41])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax2.tick_params(labelleft='off')  
    ax2.yaxis.tick_right()

    ax2.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',      # ticks along the bottom edge are off
        right='on',         # ticks along the top edge are off
        labelbottom='off',
        labelright='off'
        )


    plt.legend(loc="upper right")
    plt.xlabel('Clique Size')
    plt.ylabel('F1-score')

    plt.savefig(plot_filename)


if __name__ == "__main__":
    main()



