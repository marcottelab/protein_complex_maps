
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import argparse

import pandas as pd
import itertools as it
import scipy
from scipy.stats import pearsonr
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'grid': False})
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

from collections import Counter

import protein_complex_maps.protein_util as pu

def main():
    parser = argparse.ArgumentParser(description="Plot APMS matrix")
    parser.add_argument("--cluster_file", action="store", dest="cluster_file", required=True, 
                                            help="Filename of clusters, geneids")
    parser.add_argument("--data_filename", action="store", dest="data_filename", required=True, 
                                            help="Filename of data (format: index, experiment_id, geneid, abundance, dataset)")
    parser.add_argument("--cluster_ids", action="store", dest="cluster_ids", nargs='+', required=True, 
                                            help="Cluster ids in which to plot")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename of output plot")
    parser.add_argument("--mapping_file", action="store", dest="mapping_file", required=False, default=None,
                                            help="Filename of mapping from geneid to acc to genename")
    parser.add_argument("--map_datasets", action="store", dest="map_datasets", nargs='+', required=False, default=[], 
                                            help="List of datasets to convert experiment_ids to genenames")
    parser.add_argument("--baits_in_cluster_only", action="store_true", dest="baits_in_cluster_only", required=False, default=False, 
                                            help="Only show baits present in cluster, default=False")
    parser.add_argument("--clustering_metric", action="store", dest="clustering_metric", required=False, default='correlation', 
                                            help="Clustering metric, default=correlation")

    args = parser.parse_args()

   #kdrew: read in raw data
    data_df = pd.read_csv(args.data_filename)

    #kdrew: updating "unnamed" column name to "id"
    data_df.columns = ['id'] + [x for x in data_df.columns[1:]]
    print data_df.columns

    cfile = open(args.cluster_file,"rb")
    clusters = []
    for line in cfile.readlines():
        clusters.append(line.split())


    clusters_df = pd.DataFrame()
    for i in args.cluster_ids:
        cluster_df = data_df[data_df['geneid'].astype(str).isin(clusters[int(i)])]
        cluster_df['clusterid'] = i
        clusters_df = pd.concat([clusters_df,cluster_df])

    clusters_df['experiment_id'] = clusters_df['experiment_id'].astype(int).astype(str)


    if args.mapping_file != None:
        mapping_df = pd.read_csv(args.mapping_file,sep="\t")
        #kdrew: adjust columns and parse id sets
        mapping_df.columns = ['acc','genename_set','geneid_set']
        mapping_df['geneid_set'] = mapping_df['geneid_set'].fillna('')
        mapping_df['genename_set'] = mapping_df['genename_set'].fillna('')
        mapping_df['geneid_map'] = mapping_df['geneid_set'].apply(lambda k: k.split(';')[0])
        mapping_df['genename'] = mapping_df['genename_set'].apply(lambda k: k.split(' ')[0])


        #kdrew: merge with clusters dataframe
        clusters_df_merge = clusters_df.merge(mapping_df, how="left", left_on="geneid", right_on="geneid_map")
        #kdrew: set id to be index for merging 
        clusters_df_merge.set_index('id', inplace=True)


        for ds in args.map_datasets:
            #kdrew: for each inputted dataset, merge mapping
            clusters_df_exp_merge = clusters_df[(clusters_df.dataset == ds)].merge(mapping_df, how="left", left_on="experiment_id", right_on="geneid_map")

            #kdrew: set id to be index for updating full matrix
            clusters_df_exp_merge.set_index('id', inplace=True)
            #kdrew: update experiment_id in full matrix to be genename
            clusters_df_merge['experiment_id'].update(clusters_df_exp_merge['genename'])

        #kdrew: make pivot with genename as index
        clusters_pivot = pd.pivot_table(clusters_df_merge, values='abundance',rows=['clusterid','genename'],cols=['dataset','experiment_id'])

    else:
        clusters_pivot = pd.pivot_table(clusters_df, values='abundance',rows=['clusterid','geneid'],cols=['dataset','experiment_id'])

    sns.set_context("paper", font_scale=0.6, rc={"lines.linewidth": 2.5})

    #kdrew: create linkage matrix outside of clustermap, didn't seem to work as expected, clustering was off 
    #row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='single') for x in (clusters_pivot.values, clusters_pivot.values.T))
    #sns.clustermap(clusters_pivot.fillna(0.0),col_cluster=True,row_cluster=False,col_linkage=col_linkage)

    clusters_pivot = clusters_pivot.fillna(0.0)
    #kdrew: sort columns by sum (there has to be an easier way to do this) 
    clusters_pivotT = clusters_pivot.T
    clusters_pivotT['exp_sum'] = clusters_pivot.sum()
    clusters_pivot =  clusters_pivotT.sort('exp_sum',ascending=False).T
    clusters_pivot.drop('exp_sum', inplace=True)

    #kdrew: setup colors for rows and columns based on clusterid and dataset respectively
    #current_palette = sns.color_palette("Paired")
    current_palette = sns.color_palette(['#d6604d','#fddbc7','#d1e5f0','#4393c3','#053061'])

    if len(args.map_datasets) > 0:
        col_colors = [current_palette[args.map_datasets.index(x[0])] for x in clusters_pivot.T.index]
    else:
        col_colors = [current_palette[i] for i, x in enumerate(clusters_pivot.T.index)]

    #clusters_palette = sns.color_palette("Paired")
    clusters_palette = sns.color_palette(['#67001f','#d6604d','#fddbc7','#d1e5f0','#4393c3'])

    row_colors = [clusters_palette[-1*(1+args.cluster_ids.index(x[0]))] for x in clusters_pivot.index]

    cmap = sns.cubehelix_palette(as_cmap=True, start=2, rot=0.0, dark=0.3, light=0.95)
    #cmap = sns.cubehelix_palette(as_cmap=True)

    #kdrew: plot clustermap
    cm = sns.clustermap(clusters_pivot, col_cluster=False, cmap=cmap, col_colors=col_colors, row_colors=row_colors, row_cluster=True, cbar_kws={"label": "Percentile"}, metric=args.clustering_metric, figsize=(20, 7))

    #kdrew: setup legend for colors on rows (cluster ids) and columns (datasets)
    for i, clusterid in enumerate(args.cluster_ids):
            cm.ax_row_dendrogram.bar(0, 0, color=clusters_palette[-1*(1+i)],
                                                label=clusterid, linewidth=0)
            cm.ax_row_dendrogram.legend(loc="center left", ncol=1, fontsize=10)

    for i, dataset in enumerate(args.map_datasets):
            cm.ax_col_dendrogram.bar(0, 0, color=current_palette[i],
                                                label=dataset, linewidth=0)
            cm.ax_col_dendrogram.legend(loc="lower center", ncol=1, fontsize=10)


    #kdrew: relabel rows
    ylabels = [y.get_text().split('-')[1] for y in cm.ax_heatmap.yaxis.get_majorticklabels()]
    cm.ax_heatmap.set_yticklabels(ylabels)
    plt.setp(cm.ax_heatmap.get_yticklabels(), fontsize=10)

    #kdrew: relabel columns
    xlabels = [x.get_text().split('-')[1] for x in cm.ax_heatmap.xaxis.get_majorticklabels()]
    print Counter(xlabels)
    #kdrew: only show labels for baits that were also in clusters
    if args.baits_in_cluster_only:
        xlabels = [x if x in ylabels else ' ' for x in xlabels]
    cm.ax_heatmap.set_xticklabels(xlabels)
    plt.setp(cm.ax_heatmap.get_xticklabels(), fontsize=10)

    #kdrew: remove labels on clustermap axes
    cm.ax_heatmap.set_xlabel('')
    cm.ax_heatmap.set_ylabel('')

    if args.plot_filename is None:
        print "plot_filename is None"
        plt.show()
    else:
        plt.savefig(args.plot_filename)


if __name__ == "__main__":
	main()



