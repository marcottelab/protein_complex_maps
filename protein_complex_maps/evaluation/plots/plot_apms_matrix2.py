from __future__ import print_function

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
#rcParams.update({'grid': False})
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

from collections import Counter

import protein_complex_maps.util.protein_util as pu

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
    parser.add_argument("--no_xlabels", action="store_true", dest="no_xlabels", required=False, default=False, 
                                            help="Do not show x labels for each experiment, default=False")
    parser.add_argument("--stagger", action="store_true", dest="stagger", required=False, default=False, 
                                            help="Stagger y labels for each experiment, default=False")
    parser.add_argument("--clustering_metric", action="store", dest="clustering_metric", required=False, default='correlation', 
                                            help="Clustering metric, default=correlation")
    parser.add_argument("--plot_protein_per_exp_limit", action="store", type=int, dest="plot_protein_per_exp_limit", required=False, default=1, 
                                            help="Limit (int) of proteins per experiment to be plotted (default = 1)")

    args = parser.parse_args()

    #kdrew: read in raw percentile data
    data_df = pd.read_csv(args.data_filename,dtype={'geneid':str,'experiment_id':str})

    #kdrew: updating "unnamed" column name to "id"
    data_df.columns = ['id'] + [x for x in data_df.columns[1:]]
    print(data_df.columns)
    #print(data_df)

    cfile = open(args.cluster_file,"r")
    clusters = []
    for line in cfile.readlines():
        clusters.append(line.split())

    #print("clusters")
    #print(clusters)


    #kdrew: clusters_df is a dataframe which holds the percentile data for the user specified clusters
    clusters_df = pd.DataFrame()
    for i in args.cluster_ids:
        cluster_df = data_df[data_df['geneid'].astype(str).isin(clusters[int(i)])]
        cluster_df['clusterid'] = i
        clusters_df = pd.concat([clusters_df,cluster_df])

    print("clusters_df")
    print(clusters_df)

    #clusters_df['experiment_id'] = clusters_df['experiment_id'].astype(int).astype(str)


    if args.mapping_file != None:
        mapping_df = pd.read_csv(args.mapping_file,sep="\t")
        #kdrew: adjust columns and parse id sets
        mapping_df.columns = ['acc','genename_set','geneid_set']
        mapping_df['geneid_set'] = mapping_df['geneid_set'].fillna('')
        mapping_df['genename_set'] = mapping_df['genename_set'].fillna('')
        mapping_df['geneid_map'] = mapping_df['geneid_set'].apply(lambda k: k.split(';')[0])
        mapping_df['genename'] = mapping_df['genename_set'].apply(lambda k: k.split(' ')[0])

        print("mapping_df")
        print(mapping_df)

        #kdrew: merge with clusters dataframe
        clusters_df_merge = clusters_df.merge(mapping_df, how="left", left_on="geneid", right_on="geneid_map")
        #kdrew: set id to be index for merging 
        clusters_df_merge.set_index('id', inplace=True)

        print("clusters_df_merge")
        print(clusters_df_merge)

        df_list = []
        for ds in args.map_datasets:
            print( "ds: %s\n" % ds)
            print(clusters_df.columns)
            #print(mapping_df.columns)
            #print("clusters_df ds:")
            #print(clusters_df[(clusters_df.dataset == ds)])

            #kdrew: for each inputted dataset, merge mapping
            clusters_df_exp_merge = clusters_df[(clusters_df.dataset == ds)].merge(mapping_df, how="left", left_on="experiment_id", right_on="geneid_map")

            #print("clusters_df_exp_merge")
            #print(clusters_df_exp_merge)

            clusters_df_exp_merge = clusters_df_exp_merge[['id','genename']].drop_duplicates(subset='id')
            #kdrew: set id to be index for updating full matrix
            clusters_df_exp_merge.set_index('id', inplace=True)

            #kdrew: slice out current dataset
            clusters_df_merge_ds = clusters_df_merge[clusters_df_merge['dataset'] == ds].drop_duplicates()
            print("clusters_df_merge_ds")
            print(clusters_df_merge_ds)
            print("clusters_df_exp_merge")
            print(clusters_df_exp_merge)
            print(len(clusters_df_exp_merge.index.values))
            print(len(set([x for x in clusters_df_exp_merge.index.values])))
            #kdrew: update experiment_id in full matrix to be genename
            clusters_df_merge_ds['experiment_id'].update(clusters_df_exp_merge['genename'])
            pd.set_option('display.max_rows', len(clusters_df_exp_merge))
            #print(clusters_df_merge_ds)
            pd.reset_option('display.max_rows')
            df_list.append(clusters_df_merge_ds)

        #kdrew: recombine individual datasets
        clusters_df_merge = pd.concat(df_list)
        print(clusters_df_merge)

        #pd.set_option('display.max_rows', len(clusters_df_merge))
        #pd.set_option('display.max_columns', 500)
        #print(clusters_df_merge)
        #pd.reset_option('display.max_rows')
        #pd.reset_option('display.max_columns')
        #kdrew: make pivot with genename as index
        #clusters_pivot = pd.pivot_table(clusters_df_merge, values='abundance',rows=['clusterid','genename'],cols=['dataset','experiment_id'])
        clusters_pivot = pd.pivot_table(clusters_df_merge, values='abundance',index=['clusterid','genename'],columns=['dataset','experiment_id'])

    else:
        clusters_pivot = pd.pivot_table(clusters_df, values='abundance',rows=['clusterid','geneid'],cols=['dataset','experiment_id'])

    #sns.set_context("paper", font_scale=1.0, rc={"lines.linewidth": 2.5})

    #kdrew: create linkage matrix outside of clustermap, didn't seem to work as expected, clustering was off 
    #row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='single') for x in (clusters_pivot.values, clusters_pivot.values.T))
    #sns.clustermap(clusters_pivot.fillna(0.0),col_cluster=True,row_cluster=False,col_linkage=col_linkage)

    clusters_pivot = clusters_pivot.fillna(0.0)
    #kdrew: sort columns by sum (there has to be an easier way to do this) 
    clusters_pivotT = clusters_pivot.T
    clusters_pivotT['exp_sum'] = clusters_pivot.sum()
    clusters_pivot =  clusters_pivotT.sort_values('exp_sum',ascending=False).T
    clusters_pivot.drop('exp_sum', inplace=True)

    #kdrew: only show experiments that have a certain limit of positive entries, limit set by user
    clusters_pivot = clusters_pivot.iloc[:,((clusters_pivot > 0.0).sum() >= args.plot_protein_per_exp_limit).values]

    #kdrew: setup colors for rows and columns based on clusterid and dataset respectively
    #current_palette = sns.color_palette("Paired")
    #current_palette = sns.color_palette(['#d6604d','#fddbc7','#d1e5f0','#4393c3','#053061'])
    #current_palette = sns.color_palette(['#da534a','#58a5a7','#dbb429','#4393c3','#053061','#61055f','#010304','#c37243'])
    current_palette = sns.color_palette(['#da534a','#58a5a7','#dbb429','#c34394','#b243c3','#61055f','#010304','#c37243'])

    if len(args.map_datasets) > 0:
        col_colors = [current_palette[args.map_datasets.index(x[0])] for x in clusters_pivot.T.index]
    else:
        col_colors = [current_palette[i] for i, x in enumerate(clusters_pivot.T.index)]

    #clusters_palette = sns.color_palette("Paired")
    #clusters_palette = sns.color_palette(['#67001f','#d6604d','#fddbc7','#d1e5f0','#4393c3'])
    clusters_palette = sns.color_palette(['#d1e5f0'])

    row_colors = [clusters_palette[-1*(1+args.cluster_ids.index(x[0]))] for x in clusters_pivot.index]

    cmap = sns.cubehelix_palette(as_cmap=True, start=2, rot=0.0, dark=0.3, light=0.95)
    #cmap = sns.cubehelix_palette(as_cmap=True)

    #kdrew: plot clustermap
    print("clusters_pivot")
    #print(clusters_pivot)
    print(clusters_pivot.iloc[:,((clusters_pivot > 0.0).sum() > 1).values])

    #cm = sns.clustermap(clusters_pivot, col_cluster=False, cmap=cmap, col_colors=col_colors, row_colors=row_colors, row_cluster=True, cbar_kws={"label": "Rank Percentile"}, metric=args.clustering_metric, figsize=(20, 10))
    cm = sns.clustermap(clusters_pivot, col_cluster=True, cmap=cmap, col_colors=col_colors, row_colors=row_colors, row_cluster=True, cbar_kws={"label": "Rank Percentile"}, metric=args.clustering_metric, figsize=(20, 10))
    #cm = sns.clustermap(clusters_pivot, col_cluster=False, cmap=cmap, col_colors=col_colors, row_colors=row_colors, row_cluster=True, cbar_kws={"label": "Rank Percentile"}, metric=args.clustering_metric, figsize=(len(clusters_pivot.columns)/10,len(clusters_pivot.index))) 
    #cm = sns.clustermap(clusters_pivot, col_cluster=False, cmap=cmap, col_colors=col_colors, row_colors=row_colors, row_cluster=True, cbar_kws={"label": "Rank Percentile"}, metric=args.clustering_metric)

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
    plt.setp(cm.ax_heatmap.get_yticklabels(), fontsize=8)
    #kdrew: some reasons this doesn't work
    plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0)
    cm.ax_heatmap.yaxis.tick_left()

    print("####")
    for x in cm.ax_heatmap.xaxis.get_majorticklabels():
        print(x.get_text())
    print("####")

    #kdrew: relabel columns
    xlabels = [x.get_text().split('-')[1] for x in cm.ax_heatmap.xaxis.get_majorticklabels()]
    #print(Counter(xlabels))
    #kdrew: only show labels for baits that were also in clusters
    if args.baits_in_cluster_only:
        xlabels = [x if x in ylabels else ' ' for x in xlabels]
    if args.no_xlabels:
        xlabels = [' ' for x in xlabels]

    #kdrew: getting super hacky here
    #kdrew: for small number of experiments the columns are large so there is no overlap of labels 
    #kdrew: but requires staggering if number of experiments grows large (picking 100 was a guess)
    if len(xlabels) < 100:
        labels_in_row_threshold  = 0
    else:
        labels_in_row_threshold = 6

    #kdrew: for small clusters (<5 subunits), use small initial offset (keeps xlabels close to plot), otherwise make larger offset so labels don't overlap plot
    if len(ylabels) < 5:
        initial_offset = -0.1
    else:
        initial_offset = -1.0

    if args.stagger:
        #kdrew: annotate column labels 
        #kdrew: stagger moves the text lower
        #kdrew: empties_in_row keeps track of when to move labels back closer to start label position
        stagger = 0
        empties_in_row = 0
        labels_in_row = 0
        for i,x in enumerate(cm.ax_heatmap.xaxis.get_majorticklocs()):
            #print("stagger, empties, labels")
            #print(stagger)
            #print(empties_in_row)
            #print(labels_in_row)
            if xlabels[i] != ' ':
                empties_in_row = 0
                #cm.ax_heatmap.annotate(xlabels[i], xy=(x, 0), xytext=(x, -0.5-stagger), rotation=90, fontsize=4,)
                cm.ax_heatmap.annotate(xlabels[i], xy=(x, 0), xytext=(x, initial_offset-stagger), rotation=90, fontsize=10,)
                        #arrowprops=dict(facecolor='black', width=0.25, headwidth=0.5, shrink=0.05),)
                stagger+= 1+.1*len(xlabels[i])
                labels_in_row+=1

                #kdrew: reset labels in row
                if labels_in_row > labels_in_row_threshold:
                    stagger = 0
                    labels_in_row = 0

            else:
                labels_in_row = 0
                empties_in_row+=1
                if empties_in_row > 4:
                    stagger = 0
        cm.ax_heatmap.set_xticklabels([' ' for x in xlabels])
    else:
        cm.ax_heatmap.set_xticklabels(xlabels)
        plt.setp(cm.ax_heatmap.get_xticklabels(), fontsize=10)

    #kdrew: remove labels on clustermap axes
    cm.ax_heatmap.set_xlabel('')
    cm.ax_heatmap.set_ylabel('')

    #plt.gcf().subplots_adjust(bottom=0.45)

    if args.plot_filename is None:
        print("plot_filename is None")
        plt.show()
    else:
        #plt.savefig(args.plot_filename, dpi=300)
        cm.savefig(args.plot_filename, dpi=300)


if __name__ == "__main__":
	main()



