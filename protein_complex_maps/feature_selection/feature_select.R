#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
args <- commandArgs(trailingOnly = TRUE)

print(args)

print(args[1])

#features <- read.csv("arathtraesorysj_euNOG_corum_train_feature_selection.csv", sep=",")
features <- read.csv(args[1], sep=",")

print(head(features))


#levels(features$feature) <- order(features$feature)

f <- ggplot(features, aes(x=feature, y=score)) +
    geom_bar(stat="identity") +
    geom_text(aes(x=feature,y =score, label=rank, vjust=0)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
    background_grid(major = "y", minor = "none") +
    facet_wrap(~model)

feature_plot <- paste("featureplot_", args[1], ".png", sep="")
ggsave(filename=feature_plot, width=15, height=15, units="in")



#RandomForestRegressor
RFR_outfile = paste("features_RFR_", args[1], sep="")
RFR_features <- features %>% filter(model == "RandomForestRegressor") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(RFR_features, file=RFR_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#RandomForestClassifier
RFC_outfile = paste("features_RFC_", args[1], sep="")
RFC_features <- features %>% filter(model == "RandomForestClassifier") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(RFC_features, file=RFC_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#AdaBoostClassifier
ABC_outfile = paste("features_ABC_", args[1], sep="")
ABC_features <- features %>% filter(model == "AdaBoostClassifier") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(ABC_features, file=ABC_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#AdaBoostRegressor
ABR_outfile = paste("features_ABR_", args[1], sep="")
ABR_features <- features %>% filter(model == "AdaBoostRegressor") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(ABR_features, file=ABR_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#GradientBoostingClassifier
GBC_outfile = paste("features_GBC_", args[1], sep="")
GBC_features <- features %>% filter(model == "GradientBoostingClassifier") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(GBC_features, file=GBC_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#GradientBoostingRegressor
GBR_outfile = paste("features_GBR_", args[1], sep="")
GBR_features <- features %>% filter(model == "GradientBoostingRegressor") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(GBR_features, file=GBR_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#ExtraTreesRegressor
RFR_outfile = paste("features_RFR_", args[1], sep="")
RFR_features <- features %>% filter(model == "ExtraTreesRegressor") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(RFR_features, file=RFR_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

#ExtraTreesClassifier
RFC_outfile = paste("features_RFC_", args[1], sep="")
RFC_features <- features %>% filter(model == "ExtraTreesClassifier") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(RFC_features, file=RFC_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)


#RandomTreesEmbedding
RTE_outfile = paste("features_RTE_", args[1], sep="")
RTE_features <- features %>% filter(model == "RandomTreesEmbedding") %>% filter(rank %in% c(0,1,2,3,4,5,6,7,8,9)) %>% select(feature) 
write.table(RTE_features, file=RTE_outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)


