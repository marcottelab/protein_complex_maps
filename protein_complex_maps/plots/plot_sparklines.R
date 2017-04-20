#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(purrr))

suppressPackageStartupMessages(library(cowplot))

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-f", "--filenames", nargs='+', dest="filenames", required=TRUE,
                        help="Input elution profile")
parser$add_argument("-o", "--output_filename", action="store", required=TRUE,
                        dest="output_filename", help="Output filename for resulting plot")
parser$add_argument("-p", "--proteins", nargs='+', dest="protein_ids", required=TRUE,
                        help="Proteins to plot")
parser$add_argument("-n", "--parse_fraction_name", nargs='+', dest="parse_fraction_name", required=TRUE,
                        help="Fields to separate fraction name into, example: cell_type condition col_type fraction date, use -n fraction for default behavior")
parser$add_argument("-s", "--fraction_name_sep", dest="fraction_name_sep", required=TRUE,
                        help="Separator to separate fraction name, (would set default but R's argparse has a bug)")
parser$add_argument("-m", "--id_mapping", dest="id_mapping", required=TRUE,
                        help="Uniprot mapping file, (would set default but R's argparse has a bug)")
parser$add_argument("-g", "--id_column", dest="id_column", required=TRUE,
                        help="Column to use for protein identifier, (would set default but R's argparse has a bug)")

args <- parser$parse_args()

#filename <- "/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_helaN_1003.prot_count_uniqpeps2_FDR0010"
#output_filename <- "/project/kdrew/data/tmp/hela_tidy_trim_df.pdf"
#protein_ids <- c('ENSG00000000460','ENSG00000001167','ENSG00000002746','ENSG00000002834')

normalit<-function(m){
    (m - min(m))/(max(m)-min(m))
}

#kdrew from CMcWhite
bare_sparklines <- function(z, return_leg=FALSE){
    #z2 <- z %>% group_by(experiment) %>% mutate(speccounts = normalit(speccounts))
    #plt <- ggplot(z, aes(x=fraction, y=speccounts, group=ID, color=ID)) +
    plt <- ggplot(z, aes(x=fraction, y=speccounts, group=condition, color=condition, linetype=condition )) +
        #kdrew: set the y axis label to be protein id
        ylab(z[[args$id_column]][1]) +
        geom_line(alpha=0.9) +
        scale_linetype_manual(values=c("solid","F1")) + 
        scale_color_manual(values= c("black", "firebrick2","black","dodgerblue4","#8b104e")) 
    #kdrew: return 
    if(return_leg){
        plt <- plt + theme(legend.title=element_blank(), legend.text=element_text(size=6), legend.position="bottom" )
    }
    else{
        plt <- plt + theme(axis.text.x=element_blank(), axis.title.y = element_text(angle = 0, size=6), legend.position="none")
    }
    final_plt <- ggplotGrob(plt)
    final_plt <- gtable_remove_grobs(final_plt, c('xlab-b', 'axis-b','axis-l', 'spacer'))
    final_plt <- gtable_squash_rows(final_plt, c(1,2,3, 4, 5, 8,9, 10))
    final_plt <- gtable_squash_cols(final_plt, c(2,3,5,6))
    final_plt$vp = grid::viewport(height=0.9, width=0.9)
    return(final_plt)
}

fractionation_tidy_full_df <- NULL

for( fname in args$filenames )
{
    #kdrew: readin elution matrix
    fractionation_df = read.table(fname, header=TRUE, comment.char='!')

    #kdrew: rename protein id column
    fractionation_df$ID <- fractionation_df$X.ProtID
    fractionation_df$X.ProtID <- NULL

    #kdrew: dense matrix to sparse matrix format (tidy)
    fractionation_tidy_df <- fractionation_df %>% droplevels %>% subset(select=-TotalCount) %>% tidyr::gather(key = fraction, value = speccounts, -ID)

    print(args$parse_fraction_name)
    if(args$parse_fraction_name[1] != "fraction") {
        fractionation_tidy_df <- fractionation_tidy_df %>% separate(fraction, sep=args$fraction_name_sep, into=args$parse_fraction_name)
    }
    else{
        fractionation_tidy_df <- fractionation_tidy_df %>% mutate(condition='default')
    }


    if(is.null(fractionation_tidy_full_df)){
        fractionation_tidy_full_df <- fractionation_tidy_df
    }
    else{
        fractionation_tidy_full_df <- bind_rows(fractionation_tidy_full_df,fractionation_tidy_df)
    }
}

mapping_df <- read.table(args$id_mapping, header=TRUE, sep=',')

#kdrew: filter to just the proteins we are interested in
fractionation_tidy_trim_df <- fractionation_tidy_full_df %>% filter(ID %in% args$protein_ids)
#fractionation_tidy_trim_df_dummy <- head(fractionation_tidy_trim_df,1) #%>% filter(ID == args$protein_ids[1])
fractionation_tidy_trim_df_dummy <- fractionation_tidy_trim_df %>% filter(ID == fractionation_tidy_trim_df$ID[1])
print(tbl_df(fractionation_tidy_trim_df))
print(fractionation_tidy_trim_df$fraction)

fractionation_tidy_trim_df <- fractionation_tidy_trim_df %>% left_join(mapping_df, by=c("ID"="Entry"))

#kdrew: plot sparklines
fractionation_tidy_trim_df %>% droplevels %>% split(.$ID) %>% map(bare_sparklines) -> plotlist2
plt <- bare_sparklines(fractionation_tidy_trim_df_dummy, return_leg=TRUE)
legend <- get_legend(plt)
clusterplot <- plot_grid(plotlist = plotlist2, legend, ncol=1, align = "v")
ggsave(args$output_filename, plot=clusterplot)
#pdf(args$output_filename)
#print(clusterplot)
##kdrew: output plot to file
#dev.off()


