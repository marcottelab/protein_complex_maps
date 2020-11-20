import argparse
#from itertools import izip
#from igraph import Graph
import pandas as pd
import rpy2.robjects.packages as rpackages
utils = rpackages.importr('utils')
from rpy2.robjects.vectors import StrVector
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr

import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
 
# select a mirror for R packages
#utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R package names
#packnames = ('igraph', 'dendextend', 'dplyr')
   
#names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
#if len(names_to_install) > 0:
#    utils.install_packages(StrVector(names_to_install))

def main():

    parser = argparse.ArgumentParser(description="Cluster elutions rows based on edge scores. If weights are provided, uses 1/weight as distance")
    parser.add_argument("--input_edges", action="store", dest="input_edges", required=True, 
                                    help="Filename of edge table must be either two columns (ID1\tID2) or three columns with numeric score (ID1\tID2\tweight)(Higher score = better edge)(Default no header)")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, 
                                    help="Base name of outfile. will write both .xlsx and .csv") 
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='\t',                                                                                        help="Separator for input edge table, default=\t")
    parser.add_argument("--id_cols", action="store", nargs="+", dest="id_cols", required=False, default =["ID1", "ID2"],
                                    help="What column(s) contain identifiers. If both are in one column, provide id_sep argument, default='P1 P2'")
 
    parser.add_argument("--id_sep", action="store", dest="id_sep", required=False,
                                    help="If both identifiers are in one column with a separator")
    parser.add_argument("--weight_col", action="store", dest="weight_col", required=False,
                                    help="Column that contains column to weight network by. Higher value is stronger connection.")
    parser.add_argument("--export_excel", action = "store_true", required = False, default = False,
                                    help = "Optional flag to export an xlsx formatted sheet, may fail with encoding error")

    parser.add_argument("--header", action="store_true", dest="header", required=False, default=True,
                                    help="Flag if input has header, if no header first two columns must be identfiers") 
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
    parser.add_argument("--tree_cutoff_fractions", action="store", dest="tree_cutoff_fractions", required=False, default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1], 
                                    help="Positions for cutting the hierarchical tree (proportion of tree height), default=[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]")
    parser.add_argument("--input_elution", action="store", dest="input_elution", required=False, 
                                    help="Optional: If present, this elutions will be annotated with dendrogram order (.csv format)")
    parser.add_argument("--annotation_file", dest='annotation_file', required=False,  nargs='?', const='', 
                                    help="Optional: Path to annotations file. A tab-delimited file with protein ids in first column and annotations in second")
 
    args = parser.parse_args()


    if args.header == True:
        scores_nodups = pd.read_csv(args.input_edges, sep=args.sep)

        if len(args.id_cols) == 1:
           if not args.id_sep:
              print("If both identifiers are in one column, args.id_sep must be provided")
              return
           else:
              print(scores_nodups)
              scores_nodups[["ID1", "ID2"]] = scores_nodups[args.id_cols[0]].str.split(' ', expand=True)
              scores_nodups = scores_nodups.drop(args.id_cols, axis = 1 )


        elif len(args.id_cols) == 2:
              scores_nodups['ID1','ID2'] = scores_nodups[args.id_cols]

        else: 
             print("No more than two columns of identifiers accepted")
             return

    if args.weight_col:
            scores_nodups['weight'] = scores_nodups[args.weight_col]
            scores_nodups = scores_nodups.drop([args.weight_col], axis = 1)           
 

    else:    
         scores_nodups = pd.read_csv(args.input_edges, sep=args.sep, header=None)
         if len(scores_nodups.columns) == 2:
             scores_nodups.columns = ['ID1', 'ID2'] 
         elif  len(scores_nodups.columns) == 3:   
             scores_nodups.columns = ['ID1', 'ID2', 'weight']  
         else:
             print("Files provided with no header cannot have more than 3 columns (ID1, ID2, optional weight")
             return   

   

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
         
    if args.cutoff:
        scores_nodups = scores_nodups.head(args.cutoff)


    print(scores_nodups.head)

    #CDM: There is some problem with the python igraph dendrograms. Can't retrieve tip names or cut the dendrogram
    #graph = Graph.TupleList(scores_nodups.itertuples(index=False), weights=True, directed=False)
    #dendrogram = graph.community_walktrap(weights = 'weight')
    #dendrogram = fix_dendrogram(graph, dendrogram)
    #Gave up and imported working R code, also from igraph
    
    #Test comment out
    try:
        igraph = importr('igraph')
        dendextend = importr('dendextend')
        dplyr = importr('dplyr', on_conflict="warn")
    except Exception as e:
        print(e)
        print("Make sure igraph, dendextend, and dplyr are installed in your local R environment")
        print("Install at command line with for example:")
        print("R -e \"install.packages('igraph')\"")
        exit(1)
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
        print("Get height")
        node_attr <- get_nodes_attr(dendrogram, "height")
        print(node_attr)
        print("Get max height")

        ht <- max(node_attr) 
        print(ht)

        cut_clusters <- data.frame(ID = as.character())
        for (c in cuts){
          print(c)
          print(dim(cut_clusters))
          cut_dend <- cut_df(dendrogram, c*ht)
          print(cut_dend)
          cut_clusters <- merge(cut_clusters, cut_dend, all=TRUE)
   
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
    print(dir(pandas2ri))
    print(r_cut_clusters)
    try:
       pd_cut_clusters =  pandas2ri.ri2py_dataframe(r_cut_clusters)
    except Exception as e:
       # Method to convert from r to py was changed
       with localconverter(ro.default_converter + pandas2ri.converter):
          pd_cut_clusters = ro.conversion.rpy2py(r_cut_clusters)



       #pd_cut_clusters =  pandas2ri.rpy2py_dataframe(r_cut_clusters)
    pd_cut_clusters = pd_cut_clusters.set_index(['ID'])

   
    try:
        pd_order = pandas2ri.ri2py_dataframe(r_order)
    except Exception as e:
       with localconverter(ro.default_converter + pandas2ri.converter):
          pd_order = ro.conversion.rpy2py(r_order)


        #pd_order = pandas2ri.rpy2py_dataframe(r_order)

    pd_order = pd_order.set_index(['ID'])


    print("Combine to output table")
    outdf = pd_cut_clusters
    print(args.annotation_file)
    if args.annotation_file:
        print("Adding annotations")
        annots = pd.read_table(args.annotation_file, index_col = 0, encoding= 'unicode_escape')
        outdf = outdf.join(annots, how = 'left')


    if args.input_elution:
        print("Adding input elution data")
        elution = pd.read_csv(args.input_elution, index_col = 0)
        elution = elution.fillna(0)
        outdf = outdf.join(elution, how = "left")

    outdf = outdf.join(pd_order, how = "outer")

    outdf = outdf.sort_values(by = 'dendrogram_order')

    print("Writing csv file")
    csv_outfile = args.outfile + '.csv'
    outdf.to_csv(csv_outfile)

    if args.export_excel == True:
        print("Writing excel file, potentially will fail with encoding error")
        excel_outfile = args.outfile + '.xlsx'
        outdf.to_excel(excel_outfile)

if __name__ == "__main__":
    main()


