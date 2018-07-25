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

    parser = argparse.ArgumentParser(description="Cluster elutions rows based on edge scores. If weights are provided, uses 1/weight as distance")
    parser.add_argument("--input_edges", action="store", dest="input_edges", required=True, 
                                    help="Filename of edge table must be either two columns (ID1\tID2) or three columns with numeric score (ID1\tID2\tweight)(Higher score = better edge)(Default no header)")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, 
                                    help="Base name of outfile. will write both .xlsx and .csv") 
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='\t',
                                    help="Separator for input edge table, default=\t")
    parser.add_argument("--header", action="store_true", dest="header", required=False,
                                    help="Flag if input has header") 
    parser.add_argument("--cutoff", action="store", dest="cutoff", type = int, required=False,
                                    help="Take only first N rows")
    parser.add_argument("--threshold", action="store", dest="threshold", required=False, type = float, 
                                    help="Take only rows with score (third column) above N")
    parser.add_argument("--method", action = "store", dest = "method", required=False, default = "walktrap", choices = ["walktrap", "shortest_distance"],
                                    help="Method to calculate distance matrix, choices=walktrap, shortest_distance")

    parser.add_argument("--use_scores", action="store_true", dest="use_scores", required=False, default=True, 
                                    help="Use scores in clustering if scores present, default=True")
    parser.add_argument("--steps", action="store", dest="steps", required=False, type=int, default=4, 
                                    help="Number of steps fo walktrap to take, default=4")
    #parser.add_argument("--write_distance_matrix", action="store", dest="write_distance_matrix", required=False, default=True, 
    #                                help="Write table of pairwise distances")
    parser.add_argument("--tree_cutoff_fractions", action="store", dest="tree_cutoff_fractions", required=False, default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1], 
                                    help="Positions for cutting the hierarchical tree (proportion of tree height), default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]")
    parser.add_argument("--input_elution", action="store", dest="input_elution", required=False, 
                                    help="Optional: If present, annotate with dendrogram order information(.csv format)")
    parser.add_argument("--path_to_annotations", action='store', dest='path_to_annotations', required=False, default=None, 
                                    help="Optional: Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second")
 
    args = parser.parse_args()


    if args.header == True:
         scores_nodups = pd.read_csv(args.input_edges, sep=args.sep)

    else:    
         scores_nodups = pd.read_csv(args.input_edges, sep=args.sep, header=None)

    if len(scores_nodups.columns) == 3:
        scores_nodups.columns = ['ID1', 'ID2', 'weight']    
        print(scores_nodups.head)
        if args.threshold:
             print(args.threshold)
             scores_nodups = scores_nodups[scores_nodups['weight'] >= args.threshold]
        if args.method == "shortest_path":
             #This is necessary because the weight must be treated as a distance, where closer is better
             #Inverse to fix, following suggestion in https://toreopsahl.com/tnet/weighted-networks/shortest-paths/
             print("Inverting weights to distances for shortest path. Not necessary for walktrap, 'weights: The edge weights. Larger edge weights increase the probability that an edge is selected by the random walker. In other words, larger edge weights correspond to stronger connections.'")
 
             scores_nodups['weight'] = 1/scores_nodups['weight']

        if args.use_scores == False:
             scores_nodups['weight'] = 1
         

    elif len(scores_nodups.columns) == 2:
         scores_nodups.columns = ['ID1', 'ID2']
         scores_nodups['weight'] = 1

    else:
         print("Edge table must be either two columns (ID1\tID2) or three columns with numeric score (ID1\tID2\tweight)")   
         return

    print(scores_nodups.head)

    if args.cutoff:
        scores_nodups = scores_nodups.head(args.cutoff)

    print(scores_nodups)   
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
    get_dend <- function(df, method, steps=4){
                    #Function to get a dendrogram
                    set.seed(42)
                    graph = graph_from_data_frame(df, directed=FALSE) 
                    if(method == "walktrap"){
                        dendrogram <- as.dendrogram(cluster_walktrap(graph, steps=steps))  
                    }
                    if(method == "shortest_distance"){
                      ds <- distances(graph)
                      ds[is.infinite(ds)] <- 0
                      dendrogram  <- as.dendrogram(hclust(dist(ds)))
                    }


                    return(dendrogram)
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

    print("Calculate dendrogram")
    r_dendrogram = r_get_dend(scores_nodups, args.method, args.steps)
    print("Get row order from dendrogram")
    r_order = r_order_dend(r_dendrogram)
    print("Cut tree at multiple heights to get clusters of various stringency")
    r_cut_clusters = r_cut_dend(r_dendrogram, args.tree_cutoff_fractions)

    print("Convert back to pandas objects")
    pd_cut_clusters =  pandas2ri.ri2py_dataframe(r_cut_clusters)
    pd_cut_clusters = pd_cut_clusters.set_index(['ID'])
    pd_order = pandas2ri.ri2py_dataframe(r_order)
    pd_order = pd_order.set_index(['ID'])


    print("Combine to output table")
    outdf = pd_cut_clusters
    print(outdf)
    if args.path_to_annotations != None:
        print("Adding annotations")
        annots = pd.read_table(args.path_to_annotations, index_col = 0)
        outdf = outdf.join(annots, how = 'left')


    if args.input_elution != None:
        print("Adding input elution data")
        elution = pd.read_csv(args.input_elution, index_col = 0)
        elution = elution.fillna(0)
        outdf = outdf.join(elution, how = "left")

    outdf = outdf.join(pd_order, how = "outer")

    outdf = outdf.sort_values(by = 'dendrogram_order')

    print("Writing csv file")
    csv_outfile = args.outfile + '.csv'
    outdf.to_csv(csv_outfile)
    print("Writing excel file, might give weird encoding error")
    excel_outfile = args.outfile + '.xlsx'
    outdf.to_excel(excel_outfile)

if __name__ == "__main__":
    main()


