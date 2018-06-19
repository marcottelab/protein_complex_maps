import argparse
from itertools import izip
from igraph import Graph
import pandas as pd
import rpy2.robjects.packages as rpackages
utils = rpackages.importr('utils')
from rpy2.robjects.vectors import StrVector
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr

 
# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R package names
packnames = ('igraph', 'dendextend', 'dplyr')
   
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

def main():

    parser = argparse.ArgumentParser(description="Cluster elutions rows based on edge scores")
    parser.add_argument("--input_edges", action="store", dest="input_edges", required=True, 
                                    help="Filename of edge table must be either two columns (ID1\tID2) or three columns with numeric score (ID1\tID2\tscore)")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, 
                                    help="Base name of outfile. will write both .xlsx and .csv") 
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='\t',
                                    help="Separator for input edge table, default=\t")
    parser.add_argument("--header", action="store", dest="header", required=False, default=False,
                                    help="True if input table has a header, default=False") 
    parser.add_argument("--cutoff", action="store", dest="cutoff", type = int, required=False,
                                    help="Take only first N rows")
    parser.add_argument("--threshold", action="store", dest="threshold", required=False, type = float, 
                                    help="Take only rows with score (third column) above N")
    #parser.add_argument("--write_distance_matrix", action="store", dest="write_distance_matrix", required=False, default=True, 
    #                                help="Write table of pairwise distances")
    parser.add_argument("--tree_cutoff_fractions", action="store", dest="tree_cutoff_fractions", required=False, default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1], 
                                    help="Positions for cutting the hierarchical tree (proportion of tree height), default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]")
    parser.add_argument("--input_elution", action="store", dest="input_elution", required=False, 
                                    help="Optional: If present, annotate with dendrogram order information(.csv format)")
    parser.add_argument("--path_to_annotations", action='store', dest='path_to_annotations', required=False, default=None, 
                                    help="Optional: Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second")
 
    args = parser.parse_args()

    if args.header:

         scores_nodups = pd.read_csv(args.input_edges, sep=args.sep)
    else:

         scores_nodups = pd.read_csv(args.input_edges, sep=args.sep, header=None)

    if len(scores_nodups.columns) == 3:
         scores_nodups.columns = ['ID1', 'ID2', 'weight']    

    elif len(scores_nodups.columns) == 2:
         scores_nodups.columns = ['ID1', 'ID2']
         scores_nodups['weight'] = 1

    else:
         print("Edge table must be either two columns (ID1\tID2) or three columns with numeric score (ID1\tID2\tscore)")   
         return

    print(scores_nodups.head)

    if args.cutoff:
        scores_nodups = scores_nodups.head(args.cutoff)

    if args.threshold:
        scores_nodups = scores_nodups[scores_nodups['score'] >= args.threshold]
 
    #CDM: There is some problem with the python igraph dendrograms. Can't retrieve tip names or cut the dendrogram
    #graph = Graph.TupleList(scores_nodups.itertuples(index=False), weights=True, directed=False)
    #dendrogram = graph.community_walktrap(weights = 'weight')
    #dendrogram = fix_dendrogram(graph, dendrogram)
    #Gave up and imported working R code, also from igraph

    igraph = importr('igraph')
    dendextend = importr('dendextend')
    dplyr = importr('dplyr', on_conflict="warn")

    #Load custom R functions
    robjects.r('''
    get_dend <- function(df){
                    #Function to get a dendrogram
                    set.seed(42)
                    graph = graph_from_data_frame(df, directed=FALSE) 
                    wt_dendrogram <- as.dendrogram(cluster_walktrap(graph))  
                    return(wt_dendrogram)
        }
 

    order_dend <- function(dendrogram){
        #Function to get order from a dendrogram
        names_dendrogram <- labels(dendrogram)
        ID_order <- as.data.frame(names_dendrogram)
        ID_order$order <- row.names(ID_order)
        names(ID_order) <- c("ID" , "dendrogram_order")
        ID_order <- ID_order %>% 
                       mutate(dendrogram_order = as.numeric(dendrogram_order)) %>% 
                       arrange(dendrogram_order) 
        return(ID_order)
    }

     cut_df <- function(dendrogram, height){
        #Function to cut a dendrogram
        cd <- cutree(dendrogram, h = height) %>% as.data.frame()
        cd$ID <- row.names(cd)
        cd <- cd %>% as_tibble()
        colname <- paste("cut", as.character(height), sep = "_")
        names(cd) <- c(colname, "ID") 
        return(cd)
    }
 
    cut_dend <- function(dendrogram, cuts){
        #Function to cut the dendrogram at particular heights
        ht <- max(get_nodes_attr(dendrogram, "height"))   
        cut_clusters <- data.frame(ID = as.character())
        for (c in cuts){
        cut_clusters <- merge(cut_clusters, cut_df(dendrogram, c*ht), all=TRUE)
   
         }
       return(cut_clusters)
      }
    ''')

    r_get_dend = robjects.globalenv['get_dend']
    r_order_dend = robjects.globalenv['order_dend']
    r_cut_dend = robjects.globalenv['cut_dend']

    #Calculate dendrogram by walktrap method
    r_dendrogram = r_get_dend(scores_nodups)
    #Get row order from dendrogram
    r_order = r_order_dend(r_dendrogram)
    #Cut tree at multiple heights to get clusters of various stringency
    r_cut_clusters = r_cut_dend(r_dendrogram, args.tree_cutoff_fractions)

    #Convert back to pandas objects
    pd_cut_clusters =  pandas2ri.ri2py_dataframe(r_cut_clusters)
    pd_cut_clusters = pd_cut_clusters.set_index(['ID'])
    print(pd_cut_clusters.head)
    pd_order = pandas2ri.ri2py_dataframe(r_order)
    print(pd_order.head)
    pd_order = pd_order.set_index(['ID'])


    #Combine to output table
    outdf = pd_cut_clusters
    print(outdf)
    if args.path_to_annotations != None:
        print("Adding annotations")
        annots = pd.read_table(args.path_to_annotations, index_col = 0)
        outdf = outdf.join(annots, how = 'left')


    if args.input_elution != None:
        print("Adding input elution data")
        elution = pd.read_csv(args.input_elution, index_col = 0)
        outdf = outdf.join(elution, how = "left")

    print(outdf.head)
    outdf = outdf.join(pd_order, how = "outer")

    outdf = outdf.sort_values(by = 'dendrogram_order')

    excel_outfile = args.outfile + '.xlsx'
    csv_outfile = args.outfile + '.csv'
    print "Writing excel file"
    outdf.to_excel(excel_outfile)
    outdf.to_csv(csv_outfile)

if __name__ == "__main__":
    main()


