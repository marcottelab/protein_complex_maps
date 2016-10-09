




library(tidyr)
library(ggplot2)
library(dplyr)

ordering_fscore <- read.csv("/home/kdrew/programs/libsvm-3.20/tools/test_arathtraesorysjbraolselml_pval_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt.fscore", sep="\t", header=FALSE) 

feature_order <- read.csv("test_arathtraesorysjbraolselml_pval_euNOG_feature_order.txt", sep=",", header=FALSE)

#Get rid of ID1 ID2 columns in feature matrix
feature_order <- feature_order[3:ncol(feature_order)]

#transpose 
feature_order <- t(feature_order)

colnames(feature_order) <-c("feature_id")

df_feature_order <- as.data.frame(feature_order)
df_feature_order$feature <- 1:nrow(df_feature_order)
as.character(df_feature_order$feature) -> df_feature_order$feature




infile <- read.csv("test_arathtraesorysjbraolselml_pval_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt", sep=" ", header=FALSE)
names(infile)
test <- infile %>% gather(tmp, score_tmp, -V1)
head(test)
test$tmp <- NULL
test2 <- test %>% separate(score_tmp, c("feature", "score"), sep=":")

test2 <- test2[!is.na(test2$feature),]

test3 <- merge(test2, df_feature_order, on="feature")


for_order <- as.data.frame(sapply(ordering_fscore,gsub,pattern=":",replacement=""))


#DNW
#df_feature_order[match(for_order$V1, df_feature_order$feature),]


#levels(df_feature_order$feature) <- for_order$V1



#levels(df_feature_order$feature_id) <- df_feature_order$feature_id

#test3$feature <- factor(test3$feature)


#THE APPROACH
#Get d_feautre_order reordered the right way

levels(test3$feature) <- for_order$V1

test3$V1 <- factor(test3$V1)

#levels(test3$feature_id) <- for_order$V1

#uncomment this
#test3$V1 <- factor(test3$V1)



test3$score <- as.numeric(test3$score)


#Situation... The plots are not ordered by the fscore order. Otherwise correct

plt <- ggplot(test3, aes(score, group=V1, fill=V1)) + 
       geom_density(adjust=10, aes(alpha=0.5)) +
#       geom_label(aes(x= 0, y=10,label=feature_id))+
       facet_wrap(~feature_id)
plt <- ggplot(test3, aes(score, group=V1, fill=V1)) + 
       geom_histogram(binwidth=0.1, position="dodge") +
#       geom_label(aes(x= 0, y=10,label=feature_id))+
       facet_wrap(~feature_id)

     
print(plt)




