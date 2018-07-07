library(stringr)
library(tidyr)
library(argparse)
library(dplyr)
library(readr)  #for read_csv instead of read.csv
library(ggplot2)
library(cowplot) #for nice plot formattting
library(ggridges) #for sparklines


parser <- ArgumentParser(description='Make an on demand sparkline figure')

parser$add_argument("-f", "--filenames", nargs='+', dest="filenames", required=TRUE,
                        help="Input elution profile")
parser$add_argument("-o", "--output_filename", action="store", required=TRUE,
                        dest="output_filename", help="Output filename for resulting plot")
parser$add_argument("-p", "--proteins", nargs='+', dest="protein_ids", required=TRUE,
                        help="Proteins to plot")
parser$add_argument("-n", "--parse_fraction_name", nargs='+', dest="parse_fraction_name", required=TRUE,
                        help="Fields to separate fraction name into, example: cell_type condition col_type fraction date, use -n fraction for default behavior")
parser$add_argument("-m", "--id_mapping", dest="id_mapping", required=TRUE,
                        help="Uniprot mapping file")
parser$add_argument("-g", "--id_column", dest="id_column", default="genename", help="Broken, has to be genename. Column to use for protein identifier, default = genename")
parser$add_argument("-s", "--fraction_name_sep", dest="fraction_name_sep", default='_', help="Separator to separate fraction name, default = _")



args = parser$parse_args()

normalit<-function(m){
    (m - min(m))/(max(m)-min(m))
}


complex_plot <- function(tidy_elution){

        
     plt <- ggplot(data = tidy_elution, aes(x = as.numeric(as.character(fractionnum)), y = ID, height =  norm_speccounts, group = paste(ID, condition), color = condition )) +
        geom_ridgeline(stat = "identity", alpha = 0, scale = 0.9, size = 0.3) +
        theme_ridges() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values=c("solid","longdash")) +
        scale_y_discrete(c(0,1)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y= element_text(size = 8)) +
        theme(legend.position = "top") + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        #facetting options
        #theme(panel.spacing.x = unit(0.1, "lines")) + 
        #facet_wrap(~spec + experiment_type + tissue,  nrow=1, scales="free_x")+
        NULL 
     return(plt)
    
}
 
print(args)

full_exp_tidy <- data.frame(fraction = character(), speccounts = numeric(), ID = character())

for(fname in args$filenames){

   exp <- read_delim(fname, delim = "\t")
   #Format column names and tidy elution
   exp_tidy <- exp %>% 
              rename(ID = `#ProtID`) %>% 
              select(-TotalCount) %>% 
              gather(fraction, speccounts, -ID)


   if(args$parse_fraction_name[1] != "fraction") {
      exp_tidy <- exp_tidy %>% 
              separate(fraction, into = args$parse_fraction_name, sep = args$fraction_name_sep, remove = FALSE)

    }
   #Combine any read in files
   full_exp_tidy <- rbind(full_exp_tidy, exp_tidy)
}

print(full_exp_tidy)

#Filter full dataset, and normalize spectral counts between 0 and 1
#Also remove extra digits from the fractionnum column
filt_exp_tidy <- full_exp_tidy %>% 
                          filter(ID %in% args$protein_ids) %>% 
                          group_by(ID) %>%
                          mutate(norm_speccounts = normalit(speccounts)) %>%
                          ungroup %>% 
                          mutate(fractionnum = str_extract(fractionnum,"(\\d)+"))


#Add on uniprot annotation file
mapping_df <- read_delim(args$id_mapping,  ",")
filt_exp_tidy <- filt_exp_tidy %>% 
                      left_join(mapping_df, by = c("ID" = "Entry"))

#Make new column name for selected column
#print(args$id_column)
#print(filt_exp_tidy[,args$id_column])


#Fixed using fixed column genename
filt_exp_tidy$ID2 <- filt_exp_tidy$genename #filt_exp_tidy[,args$id_column] 

#print(filt_exp_tidy)
#Replace any missing IDs 
#print(filt_exp_tidy)
#filt_exp_tidy$ID2[is.na(filt_exp_tidy$ID2)] <- filt_exp_tidy$ID
#print(filt_exp_tidy)

#Overwrite the original ID column with selected ID column

filt_exp_tidy$ID <- filt_exp_tidy$ID2


print("starting plotting")
clusterplot <- complex_plot(filt_exp_tidy)

ggsave(args$output_filename, plot = clusterplot)


