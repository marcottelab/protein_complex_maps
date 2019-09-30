
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


    df_complete %>% mutate(ID = fct_relevel(ID, {{ ID_order }}))

    plt <- ggplot(df_complete, aes( x = fct_rev(fct_reorder(FractionID, FractionID_order)),
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

    filename <- paste0("static/complex_images/", clustid_set, "_", as.character(clustid), ".svg")
    print(filename)


    plt %>% save_plot(filename, ., device = "svg")

}




safe_sparklines <- safely(sparkline_fxn)

tissue_order <- c("green", "sprout", "dark", "nuclei", "seed")


fraction_details <- read_csv("static/data/Fraction_Details.csv")

experiment_order <- read_csv("static/data/Experiment_Order.csv")

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



all_fraction_info <- left_join(fraction_order, experiment_order, by= c("ExperimentID")) %>%
  left_join(species_attribs, by = "spec")



virNOG_elut_norm_annot <- virNOG_elut_norm  %>%
  left_join(all_fraction_info, by = c("FractionID", "ExperimentID", "spec"))

complex_info <- read_csv("clustid_key.csv")

complex_info   %>%
  mutate(splitter = paste0(clustid_set, clustid)) %>%
  split(.$splitter) %>%
  map(~safe_sparklines(virNOG_elut_norm_annot, all_fraction_info, ., header = TRUE))

