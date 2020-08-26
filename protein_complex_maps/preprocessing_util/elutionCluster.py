import os.path
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

def main():
    
    parser = argparse.ArgumentParser(description="Cluster a single or multiple joined elution files")
    parser.add_argument("--input_elutions", action='store', dest='input_elutions', required=True, nargs='+', help="Input csv elution files")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, help="Name of outfile. Must end with .xlsx")
    parser.add_argument("--hclust_method", nargs="+", action="store",dest="hclust_method",required=False,default='average',choices=["single","complete","average","weighted","centroid","median","ward"],help="Method for hierarchical clustering of rows")
    parser.add_argument("--hclust_metric", nargs="+", action="store",dest="hclust_metric",required=False,default="euclidean",choices=["euclidean","canberra","braycurtis","pearson","spearman"],help="Distance or correlation metric used for clustering")
    parser.add_argument("--path_to_annotations", action='store', dest='path_to_annotations', required=False, default=None, help="Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second")
    parser.add_argument("--threshold", dest="threshold", required=False, default=False, type=int, help="Drop proteins with fewer than this many PSMs across experiments")
    parser.add_argument("--normalize", dest="normalize", nargs="+", required=False, default=False, choices=["row_sum","row_max","column"], help="Normalize data. Can be one or several of: 'row_max','row_sum','column'")
    parser.add_argument("--elution_separator", action='store', dest='sep', required=False, default="\t", help="Separator used for input elution profiles, default='\t'")
    
    args = parser.parse_args()
    
    assert args.outfile.split(".")[-1] == "xlsx", "Outfile must be .xlsx - sorry"
    
    if args.path_to_annotations != None:
        assert os.path.isfile(args.path_to_annotations), "File {} not found".format(args.path_to_annotations)
        
    
    
    ## Read in and join filesd
    joined_elutions = pd.DataFrame()
    psm_sum_df = pd.DataFrame()
    for f in args.input_elutions:
        print "Reading in {}".format(f)
        df = pd.read_csv(f,index_col=0, sep=args.sep)
        if "TotalCount" in df.columns:
            df.drop("TotalCount",inplace=True,axis=1)
        
        rowSum = pd.DataFrame( df.sum(axis=1) )
        rowSum.columns = [f.split(".")[0]]
        psm_sum_df = psm_sum_df.join(rowSum, how='outer')
        
        ## Normalize, if flagged
        if args.normalize:
            print "Normalize by {}".format(", ".join(args.normalize))
            if 'column' in args.normalize:
                df = df/df.sum()
            if 'row_sum' in args.normalize:
                if "row_max" in args.normalize:
                    print "Warning: It's kinda weird to normalize both row_sum and row_max"
                df = df.div(df.sum(axis=1), axis=0)
            if 'row_max' in args.normalize:
                df = df.div(df.max(axis=1), axis=0)
            
        joined_elutions = joined_elutions.join(df, how='outer')

    joined_elutions.fillna(0,inplace=True)
    column_order = list(joined_elutions.columns) # keep track of columns for re-ordering below
            
    ## Will be added to df later
    colSums = pd.Series(joined_elutions.sum(),name="colSum")
    rowSums = psm_sum_df.sum(axis=1)
    
   # Thresholding
    if args.threshold:
        print "Removing rows with fewer than {} PSMs".format(args.threshold)
        ## Add a temporary column for thresholding, then drop it so we don't cluster on it
        joined_elutions["tmp_rowSum"] = rowSums
        before = len(joined_elutions)
        joined_elutions = joined_elutions[joined_elutions.tmp_rowSum >= args.threshold]
        after = len(joined_elutions)
        print "Removed {} rows".format(before - after)
        rowSums = joined_elutions.pop("tmp_rowSum")
            
    labels_index = pd.Series(joined_elutions.index,index=range(len(joined_elutions))) # Create indexer to manage scipy output
    row_orderings = [] # DataFrames indexed on ids that hold the clustered orders.
    clustered_names = [] # Used during column re-ordering
    
    for method in args.hclust_method:
        for metric in args.hclust_metric:
            print "Clustering method={}, metric={}".format(method, metric)

            if metric == "pearson":
                corr_mat = np.around( np.corrcoef(joined_elutions), decimals = 4)
                dist_mat = (1. - corr_mat) / 2. # BJL: Invert correlations and move to 0-1 axis, i.e. convert them to distances
                distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
                link = hclust.linkage(distArray, method=method)
            elif metric == "spearman":
                corr_mat = np.around( stats.spearmanr(joined_elutions.T)[0], decimals = 4 )
                dist_mat = (1. - corr_mat) / 2.
                distArray = dist_mat[ np.triu_indices(dist_mat.shape[0],1) ]
                link = hclust.linkage(distArray, method=method)
            elif metric in ("euclidean","canberra","braycurtis"):
                link = hclust.linkage(joined_elutions, method=method, metric=metric)
            else:
                raise Exception("Clustering metric unknown: {}".format(metric))# shouldn't be necessary b/c arparse should catch it with "choices", but whatever
            
            leaves_list = hclust.leaves_list(link)
            ordered_index = labels_index.iloc[leaves_list] # reorder index from hier clustering
        
            cluster_name = "_".join(["order",metric,method])
            clustered_names.append(cluster_name)
            row_order = pd.DataFrame({cluster_name: range(len(ordered_index))}, index=ordered_index)
            assert len(row_order) == len(joined_elutions)
            row_orderings.append(row_order)


    ## Join clustered row orders
    for order in row_orderings:
        joined_elutions = joined_elutions.join(order,how='left')
    
    joined_elutions["rowSum"] = rowSums

    column_order = ['rowSum'] + clustered_names + column_order
    
    #kdrew: create empty annotation dataframe in case no annotations are read in
    if args.path_to_annotations != None:
        print "Adding annotations"
        annots = pd.read_table(args.path_to_annotations,index_col=0)
        #print "length of annotation df: {}".format(len(annots)) # Debugging
        joined_elutions = annots.join(joined_elutions,how='right')
        column_order = list(annots.columns) + column_order
    
    # reorder columns
    reordered_df = joined_elutions[column_order]
    
    # add column sum row
    outdf = reordered_df.append(colSums) # BJL: Must add annots, then row sum, the cluster order, then col sums -- in that order
    
    print "Writing excel file"
    outdf.to_excel(args.outfile)
    
if __name__ == '__main__':
    main()
