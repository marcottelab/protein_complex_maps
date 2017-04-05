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

conversion <- read.csv(args[2], sep=",", row.names=NULL)
awk_commands <- paste("awk_selections_", args[1], sep="")

print(head(features))
print(head(conversion))

names(conversion) <- c("num", "feature")

models <- c("RandomForestClassifier", "AdaBoostClassifier", "ExtraTreesClassifier", "RandomTreesEmbedding")


features <- features %>% filter(model %in% models)


#levels(features$feature) <- order(features$feature)

f <- ggplot(features, aes(x=feature, y=score)) +
    geom_bar(stat="identity") +
    geom_text(aes(x=feature,y =score, label=rank, vjust=0)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
    background_grid(major = "y", minor = "none") +
    facet_wrap(~model, ncol=1)


feature_plot <- paste("featureplot_", args[1], ".png", sep="")
ggsave(filename=feature_plot, width=20, height=20, units="in")

#this is for awk indexing

features <- merge(features, conversion, by="feature")

print(head(features))


feature_lists <- function(features, modelname){

    #score above zero
    outfile0 <- paste("features", args[1], modelname, "zerothreshold", sep="_")
    features0 <- features %>% filter(model == modelname) %>% filter(score > 0) %>% select(feature) 
    write.table(features0, file=outfile0, row.names=FALSE, col.names=FALSE, quote=FALSE)

    outfile0num <- paste("features", args[1], modelname, "zerothresholdnum", sep="_")
    features0num <- features %>% filter(model == modelname) %>% filter(score > 0) %>% select(num) 
    print(features0num)
    write.table(features0num, file=outfile0num, row.names=FALSE, col.names=FALSE, quote=FALSE)
    features0num_inc <- as.numeric(features0num$num) + 1

    exp0 <- paste(sort(features0num_inc),collapse = ",$")
    print(exp0)
    awk_exp0 <-  paste("awk -F' ' '{print $1,$", exp0, "}' > ", modelname, "_zerothresh", sep="")
    print(awk_exp0)
    write(awk_exp0, file=awk_commands)




    #top 10% of features. Only scores above zero

    outfile10 <- paste("features", args[1], modelname, "top24threshold", sep="_")
    features10 <- features %>% filter(model == modelname) %>% filter(rank %in% seq(1, 24, by=1)) %>% select(feature) 
    write.table(features10, file=outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)

    outfile10num <- paste("features", args[1], modelname, "top24thresholdnum", sep="_")
    features10num <- features %>% filter(model == modelname) %>% filter(rank %in% seq(1, 24, by=1)) %>% select(num) 
    write.table(features10num, file=outfile10num, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    features10num_inc <- as.numeric(features10num$num) + 1


    exp10 <- paste(sort(features10num_inc),collapse = ",$")
    print(exp10)
    awk_exp10 <-  paste("awk -F' ' '{print $1,$", exp10, "}' > ", modelname, "_top24", sep="")
    print(awk_exp10)
    write(awk_exp10, file=awk_commands, append=TRUE)

}

for(m in models){

feature_lists(features, m)

}


