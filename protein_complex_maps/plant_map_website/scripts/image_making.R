
library(tidyverse)
library(cowplot)
library(ggridges)
library(svglite)

#ggplot2::margin because  randomForest::margin keeps interfering
theme_cowplot_consistent_text <- function (font_size = 8) {

  theme_cowplot() %+replace%
       theme(strip.text = element_text(size = font_size),
             axis.text = element_text(colour = "black", size = font_size),
             plot.title = element_text(size = font_size),
             axis.title = element_text(size = font_size),
             legend.text = element_text(size = font_size),
             legend.title = element_text(size = font_size),
             legend.key.size = unit(0.5, "lines"),
             axis.title.x = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), vjust = 1),
             axis.text.x = element_text(margin = ggplot2::margin(t = 1, r = 0, b = 0, l = 0), vjust = 1),
             axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), angle = 90, vjust = 1),
             axis.text.y = element_text(margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0), hjust = 1)
            )

}
theme_set(theme_cowplot_consistent_text(font_size = 8))


theme_sparkline <- function () {
    theme_cowplot_consistent_text() %+replace%
        theme(
           axis.text.x = element_blank(),
           axis.text.y = element_text(vjust = -1),
           axis.ticks.x = element_blank(),
           axis.line.x = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),

           panel.spacing = unit(0.1, "lines"),
           panel.background = element_rect(fill = "white"),
           plot.background=element_rect(fill="grey90"),
           ##strip.text = element_text(angle = 45, hjust = 0),
           #strip.text = element_blank(),
           strip.background = element_rect(fill = NA, color = "black"),
           axis.line.y = element_blank(),
           plot.margin = unit(c(0, 0, 0, 0), "cm")#,
           #panel.border = element_rect(color = "grey85", fill = NA, linetype = 1, size =0.5)#,
            # axis.line.x.top = element_line(color = "red", size = 0.5, linetype = 1)
          )
}




# For making a sparkline plot, given a data frame, the fractions to plot, and proteins to plot, and optional defined order, and whether to use the facet headers
sparkline_fxn <- function(df, fracinfo, sel, specchoice = NULL, order = NULL, header = TRUE){

    sel <- sel %>% arrange(order)
    selection <- sel %>% pull(OrthogroupID)
    #order <- sel %>% pull(order)
    print(selection)
    print(order)
    df_complete <- df %>% filter(ID %in% selection) %>%

            bind_rows(fracinfo ) %>% #adding on so we have all the weekdays
        complete(nesting(FractionID, FractionID_order, ExperimentID, ExperimentID_order, tissue, spec, experiment_name, spec_order),ID, fill = list(expnorm_PeptideCount = 0)) %>% filter(!is.na(ID))

    if(!is.null(specchoice)){
       df_complete <- df_complete %>% filter(spec %in% specchoice)
    }

    if(!is.null(order)){
      ID_order <- order


    }else{
      ID_order <- selection

    }


    #ID_order_quo <- enquo(ID_order)
    df_complete %>% mutate(ID = fct_relevel(ID, {{ ID_order }}))

    plt <- ggplot(df_complete, aes( x = fct_reorder(FractionID, FractionID_order),
                                    y = fct_relevel(ID, {{ ID_order }}),
                                    height = expnorm_PeptideCount,
                                    group = fct_relevel(ID, {{ ID_order }}))) +
        geom_ridgeline(fill= NA, scale = 0.8) +
        theme_sparkline() +
        theme(strip.text = element_text(size = "6", angle = 90,  vjust = 0.5, hjust=0),
              panel.spacing = unit(0.1, "lines"),
              strip.background = element_rect(fill = "white", color = "white"),
              axis.title.y = element_blank()) +

        scale_y_discrete(expand = c(0,0.02)) + #This gets rid of the bottom grey space to help with alignment, add little bit of space to stop bottom line from being cut off
        facet_grid(~fct_reorder(experiment_name, ExperimentID_order), scales = "free_x", space = "free_x") +

      NULL

    if(header == FALSE){

      plt <- plt + theme(strip.text = element_blank())

    }

    clustid <- sel %>% pull(clustid) %>% head(1)
    clustid_set <- sel %>% pull(clustid_set) %>% head(1)

    filename <- paste0("static/complexes/", clustid_set, "_", as.character(clustid), ".svg")
    print(filename)
    plt %>% ggsave(filename, ., device = "svg", width = 10, units = "in")
    #return(plt)

}




safe_sparklines <- safely(sparkline_fxn)

#map(data, ~ safe_function(.x))

#arath_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/arath_hmmer_virNOG.mapping",  "arath", TRUE)
#braol_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/braol_hmmer_virNOG.mapping", "braol", TRUE)
#chlre_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/chlre_hmmer_virNOG.mapping", "chlre", TRUE)
#sollc_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/sollc_hmmer_virNOG.mapping", "sollc", TRUE)
#soybn_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/soybn_hmmer_virNOG.mapping", "soybn", TRUE)
#cocnu_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/cocnu_hmmer_virNOG.mapping",  "cocnu")
#orysj_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/orysj_hmmer_virNOG.mapping",  "orysj", TRUE)
#chqui_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/chqui_hmmer_virNOG.mapping",  "chqui")
#wheat_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/wheat_hmmer_virNOG.mapping",  "wheat", TRUE)
#selml_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/selml_hmmer_virNOG.mapping",  "selml", TRUE)
#cerri_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/cerri_hmmer_virNOG.mapping",  "cerri")
#cansa_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/cansa_hmmer_virNOG.mapping",  "cansa")
#maize_virNOG_orthology <- read_orthology("protein_identification/eggnog_mapping/maize_hmmer_virNOG.mapping",  "maize")
#
#plant_virNOG_orthology <- rbind(arath_virNOG_orthology, braol_virNOG_orthology, soybn_virNOG_orthology, cansa_virNOG_orthology, chqui_virNOG_orthology, sollc_virNOG_orthology, cocnu_virNOG_orthology,  wheat_virNOG_orthology, orysj_virNOG_orthology, selml_virNOG_orthology, cerri_virNOG_orthology, chlre_virNOG_orthology, maize_virNOG_orthology)
#
#
##Downloaded from uniprot, Entry Mass Length table

tissue_order <- c("green", "sprout", "dark", "nuclei", "seed")


fraction_details <- read_csv("static/data/Fraction_Details.csv")

# NEED to fix DSSO order when hopper is back
experiment_order <- read_csv("static/data/Experiment_Order.csv")

#experiment_translations <- read_csv("accessory_files/Experiment_translations.csv")

species_attribs <- read_csv("static/data/Species_attribs.csv")

virNOG_elut_norm <- read_csv("static/data/virNOG_elut_norm.csv") %>% filter(spec != "euggr")

#Fix ordering issues with certain poorly labeled fractions
virNOG_elut_norm <- virNOG_elut_norm %>%
      mutate(FractionID = case_when(FractionID == "wheatgerm_SEC2_8_1a_20150710" ~ "wheatgerm_SEC2_08_1a_20150710",
                                    FractionID == "wheatgerm_SEC2_4_1a_20150710" ~ "wheatgerm_SEC2_04_1a_20150710",
                                     FractionID == "wheatgerm_SEC2_6_1a_20150709" ~ "wheatgerm_SEC2_06_1a_20150709",
                                     TRUE ~ FractionID
                                    ))

fraction_order <- virNOG_elut_norm %>% select(ExperimentID, FractionID) %>%
     unique %>%
      arrange(FractionID) %>%
      group_by(ExperimentID) %>%
      mutate(FractionID_order = rank(FractionID)) %>%
      ungroup


#What is this line for?
experiment_order  <- experiment_order %>% select(ExperimentID, ExperimentID_order, tissue, spec, experiment_name, experiment_type, experiment_type_short, tissue_short)

all_fraction_info <- left_join(fraction_order, experiment_order, by= c("ExperimentID")) %>%
  left_join(species_attribs, by = "spec")



virNOG_elut_norm_annot <- virNOG_elut_norm  %>%
  left_join(all_fraction_info, by = c("FractionID", "ExperimentID", "spec"))
# Should not change row numbers

complex_info <- read_csv("clustid_key.csv")

#Will go faster on hopper
complex_info   %>%
   #filter(clustid == 26) %>%
#  filter(clustid > 88) %>%
  mutate(splitter = paste0(clustid_set, clustid)) %>%
  split(.$splitter) %>%
  map(~safe_sparklines(virNOG_elut_norm_annot, all_fraction_info, ., header = TRUE))

#sparkline_fxn(virNOG_elut_norm_annot, all_fraction_info, sel, header = TRUE)

#d_and_la_plot %>% ggsave("sparkline_d_and_la_plot.png", .,  device = "png", height = 2, width = 8, units = "in")

# Reduced one for the main text
#virNOG_elut_norm_ordered_dandlafilt <- virNOG_elut_norm_annot %>% filter(experiment_name %in% c("wheat_germ4", "cansa_seed1", "cocnu_nuc1", "braol_nuc1",  "maize_leaf1", "soybn_sprout2"))
#all_fraction_info_dandlafilt <- all_fraction_info %>% filter(experiment_name %in% c("wheat_germ4", "cansa_seed1","cocnu_nuc1", "braol_nuc1",  "maize_leaf1", "soybn_sprout2"))

d_and_la_filt_plot <- sparkline_fxn(virNOG_elut_norm_ordered_dandlafilt, all_fraction_info_dandlafilt, domino_and_la, order = c("LA1", "DOMINO1"))+ theme(axis.title.y = element_blank())



##For species with uniprot proteomes
#read_orthology <- function(filename, species, split = FALSE, order = NULL){
#    orthology <- read_delim(filename, delim = "\t")
#    if(split == TRUE){
#        orthology <- orthology %>% separate(ProteinID, into = c("tmp", "ProteinID", "tmp2"), sep = "[|]") %>% select(ID, ProteinID)
#    }
#orthology$spec <- species
 #return(orthology)
#}
#
##For species that aren't uniprot
#read_orthology2 <- function(filename, species, split = FALSE){
#    orthology <- read_delim(filename, delim = "\t")
#    if(split == TRUE){
#        orthology <- orthology %>% separate(ProteinID, into = c("tmp", "tmp2", "ProteinID"), sep = "[|]") %>% select(ID, ProteinID)
#    }
#orthology$spec <- species
#return(orthology)
#}

#What is the difference between "_ordered" and "_annot"?
#virNOG_elut_norm_ordered <- virNOG_elut_norm_annot  %>% select(-fracnorm_PeptideCount, -fracexpnorm_PeptideCount, -Total_PeptideCount, -abundance_ppm, -expnorm_ppm, -filename)


#Big clustergram
#wide <- virNOG_elut_norm_annot %>%
#    arrange(spec_order, ExperimentID_order) %>%
#    select(FractionID, ID, expnorm_PeptideCount) %>%
#    mutate(FractionID = fct_inorder(FractionID)) %>%
#    spread(FractionID, expnorm_PeptideCount, fill = 0 ) %>%
#    filter(grepl("ENOG4", ID)) %>%
#    filter(ID %in% obs_groups$ID)
#
#wide_tidy <- wide_ordered %>% gather(FractionID, value, -ID)
#alldata_plot <-  wide_tidy %>%
#     ggplot(aes( x = fct_inorder(FractionID), y = fct_rev(ID) , fill = value)) +
#     geom_tile() +
#
#     scale_fill_gradient(low = "#ffff66", high = "blue") +
#     #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust  =1))
#     theme_nothing()
#
#
#
## Known complexes
#tidy_postclust_known_filt <- tidy_postclust %>%
#
#      mutate(EXP = case_when(experiment_name == "maize_leaf1" ~ "Maize leaf\nSize exclusion fractions \u2192", #unicode arrow right
#                             experiment_name == "braol_nuc1" ~ "Broccoli nuclei\nIon exchange fractions \u2192"
#                             )) %>%
#
#      mutate(set = case_when(
#         ID %in% tset ~ "TSET\ncomplex",
#         ID %in% eif2b ~ "EIF2B\ncomplex",
#         ID %in% proteasome ~ "20S\nproteasome\ncomplex",
#         ID %in% mccase ~ "MCCase\ncomplex",
#        # ID %in% ahl ~ "AHL",
#         #ID %in% ef1 ~ "EF1",
#         #ID %in% nups ~ "NUP",
#         #ID %in% fy_cpsf ~ "FY/CPSF",
#         #ID %in% cyp1822 ~ "CYP18/22",
#         ID %in% etf ~ "ETF",
#         #ID %in% sec1331 ~ "SEC13/31",
#         #ID %in% carb ~ "Carbamoyl synthase",
#         ID %in% prefoldin ~ "prefoldin\ncomplex"
#
#
#         #ID %in% cct ~ "CCT"
#
#      )) %>%
#    filter(!is.na(set))
#
## Known complexes
## Get example order in clustered heatmaps for fig 1d reasons
#
#tidy_postclust_known_filt %>%
#     left_join(df_order_pearson_average, by = "ID") %>%
#     select(ID, simplename, order_pearson_average) %>% unique %>% View
##Proteasome is at 16,000/26,114 = 61% down
##ETFA 20374 = 78% down
#
##PFD6 19208 =73.5%  down
#
#
## amazing https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
#
#calloutdata_known_plot <- tidy_postclust_known_filt  %>% filter(experiment_name %in% c( "braol_nuc1",  "maize_leaf1")) %>%
#      arrange(ExperimentID_order) %>%
#      mutate(experiment_name = fct_inorder(experiment_name)) %>%
#     #ggplot(aes( x = fct_inorder(FractionID), y = fct_rev(simplename) , fill = set, alpha = value)) +
#    ggplot(aes( x = fct_inorder(FractionID), y = fct_rev(simplename) , fill = value)) +
#     geom_tile(color= "grey50") +
#     theme(axis.text.x = element_blank(),
#           legend.position = "none",
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           strip.text.x = element_text(hjust = 0, vjust = -1),
#           strip.background = element_rect(fill = NA),
#           axis.ticks.x = element_blank(),
#           axis.text.y = element_text(size = 7),
#           panel.spacing.x = unit(0.5, "lines"),
#           panel.spacing.y = unit(0.0001, "lines"),
#          strip.text.y = element_text(angle = 0, hjust = 0),
#          plot.margin = margin(0,0,0,0),
#          strip.placement = "outside") +
#     facet_grid(set~EXP, scales = "free", space = "free_y") +
#     scale_fill_gradient(low = "#ffff66", high = "blue") + #lighter #ffff7f, darker #ffff66, close to yellow #ffff4c
#     # scale_fill_viridis(direction =  -1) +
#    scale_x_discrete(position = "top") +
#     #panel_border() +
#      #scale_fill_manual(values = brewer.pal(name="Set1", n=9), guide="none") +
#      #scale_alpha_continuous(range=c(0, 1)) +
#    NULL
#
#
#diffclust_annot <- read_csv("training/unscaled_top100/walktraps/annotation_walktrapfdr10_3step.csv") # Loaded in twice, do once at top of script
#clusters_order <- diffclust_annot %>% select(ID) %>% rownames_to_column(var = "order")
#
#
#wide_postclust <- wide %>% inner_join(clusters_order, by = "ID")
#
#
#tidy_postclust <- wide_postclust %>%
#   gather(FractionID, value, -ID, -order) %>%
#   left_join(custom_genenames, by = "ID") %>%
#   left_join(all_fraction_info, by = "FractionID") %>%
#    arrange(order) %>%
#    mutate(ID = fct_inorder(ID), simplename = fct_inorder(simplename))
#
#
#tidy_postclust_multimer_filt <- tidy_postclust %>%
#      mutate(EXP = case_when(experiment_name == "soybn_sprout1" ~ "Soy sprout\nSize exclusion fractions \u2192", #unicode arrow right
#                             experiment_name == "maize_leaf1" ~ "Maize leaf\nSize exclusion fractions \u2192"#, #unicode arrow right
#                             #experiment_name == "arath_sprout1" ~ "arath_sprout1 \u2192", #WWC
#                             #experiment_name == "braol_nuc1" ~ " braol_nuc1\u2192", # WWC
#                             #experiment_name == "wheat_germ2" ~ "wheat_germ2 \u2192" # SEC
#
#                             )) %>%
#
#      mutate(set = case_when(
#         ID %in% csp ~ "CSP41", # dimer of trimers BAA - AAB
#         ID %in% unchar_tri ~ "Uncharacterized", # Multimer
#        # ID %in% hsp70 ~ "HSP70",
#
#         #ID %in% oxal ~ "oxal",
#         #ID %in% cleanunchar ~ "cleanunchar",
#         ID %in% apospory_unchar ~ "apospory_unchar", # Multimer
#         ID %in% glycleav_unchar ~ "clycleave_unchar" # Multimer signal
#
#
#         #ID %in% cct ~ "CCT"
#
#      )) %>%
#    filter(!is.na(set))
##
#
## amazing https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
#
#calloutdata_multimer_plot <- tidy_postclust_multimer_filt  %>%
#    #filter(experiment_name %in% c( "soybn_sprout1",  "maize_leaf1")) %>%
#    filter(!is.na(EXP)) %>%
#
#      arrange(ExperimentID_order) %>%
#      mutate(experiment_name = fct_inorder(experiment_name)) %>%
#    ggplot(aes( x = fct_inorder(FractionID), y = fct_rev(simplename) , fill = value)) +
#     geom_tile(color= "grey50") +
#     theme(axis.text.x = element_blank(),
#           legend.position = "none",
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           strip.text = element_text(hjust = 0, vjust = 0.9),
#           strip.background = element_rect(fill = NA),
#           axis.ticks.x = element_blank(),
#           #axis.line.x = element_blank(),
#           panel.spacing.x = unit(0.5, "lines"),
#           panel.spacing.y = unit(0.0001, "lines"),
#           strip.text.y = element_text(angle = 0, hjust = 0),
#          strip.placement = "outside") +
#     facet_grid(set~EXP, scales = "free", space = "free_y") +
#     scale_fill_gradient(low = "yellow", high = "blue") +
#     # scale_fill_viridis(direction =  -1) +
#     scale_x_discrete(position = "top") +
#     #panel_border() +
#      #scale_fill_manual(values = brewer.pal(name="Set1", n=9), guide="none") +
#      #scale_alpha_continuous(range=c(0, 1)) +
#    NULL
#calloutdata_multimer_plot
## Size calibration shows that these aren't single large clusters
#
#
#
#
#tidy_postclust_metabol_filt <- tidy_postclust %>%
#      mutate(EXP = case_when(experiment_name == "soybn_sprout1" ~ "Soy sprout\nSize exclusion fractions \u2192", #unicode arrow right
#                             experiment_name == "maize_leaf1" ~ "Maize leaf\nSize exclusion fractions \u2192" #unicode arrow right
#                             )) %>%
#      left_join(clusters) %>%
#      filter(cut_1085.4 == 4) %>%
#      mutate(set = cut_723.6) %>%
#      filter(! set %in% single_member_clusters) %>%
#      filter(set %in% c(17,45,61,71,80)) %>%
#
#    filter(!is.na(set))
#
#
## amazing https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
#
#calloutdata_metabol_plot <- tidy_postclust_metabol_filt  %>%
#    filter(experiment_name %in% c( "soybn_sprout1",  "maize_leaf1")) %>%
#
#      arrange(ExperimentID_order) %>%
#      mutate(experiment_name = fct_inorder(experiment_name)) %>%
#    ggplot(aes( x = fct_inorder(FractionID), y = fct_rev(simplename), fill = value)) +
#     geom_tile(color= "grey50") +
#     theme(axis.text.x = element_blank(),
#           legend.position = "none",
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           strip.text = element_text(hjust = 0, vjust = 0.9),
#           strip.background = element_rect(fill = NA),
#           axis.ticks.x = element_blank(),
#           #axis.line.x = element_blank(),
#           panel.spacing.x = unit(0.5, "lines"),
#           panel.spacing.y = unit(0.0001, "lines"),
#           strip.text.y = element_text(angle = 0, hjust = 0),
#          strip.placement = "outside") +
#     facet_grid(set~EXP, scales = "free", space = "free_y") +
#     scale_fill_gradient(low = "yellow", high = "blue") +
#
#     scale_x_discrete(position = "top") +
#
#    NULL
#
#
#
#Domino1 and La1


