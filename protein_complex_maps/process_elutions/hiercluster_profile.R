#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(purrr))

suppressPackageStartupMessages(library(cowplot))

suppressPackageStartupMessages(library(argparse))

# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input_filename", action="store", required=TRUE,
                        dest="input_filename", help="Input filename for elution_profile. One ID column and the rest numerical")
parser$add_argument("-o", "--output_filename", action="store", required=TRUE,
                        dest="output_filename", help="output filename for cluster order")
parser$add_argument("-c", "--cor_method", action="store", choices = c("pearson", "spearman", "kendall"),
                        dest="cor_method", default = "pearson", help="Correlation type for the distance function")

parser$add_argument("-hc", "--hclust_method", action="store", choices = c("average", "complete", "single", "ward.D", "ward.D2", "mcquitty", "median", "centroid" ),
                        dest="hclust_method", default = "average", help="Hclust type, default = average")
parser$add_argument("-id", "--id_column", dest="id_column", default="ID", help="Column to use for row identifier, default = ID. Everyother column should be numeric")
parser$add_argument("-s", "--sep", action="store",
                        dest="sep", default=',', help="elution profile separator. default = , ")


add_an_order <- function(df_cluster_on, df_get_order, id_col, name, cormethod, hclustmethod){
  #Takes:A dataframe to cluster on
  #      A dataframe to reorder by the clustering. For ex, cluster on cilia exps, reorder whole df
  #      The name of the ID column, The name of the new ordering columns
  #      The methods for correlation and hclust  

  hr <- hclust(as.dist(1-cor(t(df_cluster_on), method=cormethod)), method=hclustmethod)
  df_clust <- df_get_order[rev(hr$labels[hr$order]),] #reorder by clustering order

  df_clust[,id_col] <- row.names(df_clust) #remake ID columns
  rownames(df_clust) <- 1:nrow(df_clust) #Number rows 
  df_clust[,name] <- row.names(df_clust)  #The row numbers give clustered order
  df_order <- df_clust %>% select_(id_col, name)
  return(df_order)
}
args <- parser$parse_args()
print(args)

set.seed(42)

df_elut <- read.table(args$input_filename, header=TRUE, sep=args$sep)
row.names(df_elut) <- df_elut[,args$id_column]
df_elut[,args$id_column] <- NULL

df_elut[is.na(df_elut)] <- 0 

finalcolname <- paste("order", args$cor_method, args$hclust_method, sep="_")


order_col <- add_an_order(df_elut, df_elut, args$id_column, finalcolname, args$cor_method, args$hclust_method)

#outputname = paste(finalcolname, args$input_filename, sep="_")
order_col %>% write.csv(args$output_filename, quote=FALSE, row.names=FALSE) 




