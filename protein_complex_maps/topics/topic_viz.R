
library(ggplot2)
library(cowplot)


raw_data = read.csv("topic_performance.txt", header=FALSE)
f <- raw_data

 
levels(df$condition) <- c("n1000_topics_plants_euNOG_concat.csv", "n3000_topics_plants_euNOG_concat.csv", "n5000_topics_plants_euNOG_concat.csv", "n8000_topics_plants_euNOG_concat.csv", "n10000_topics_plants_euNOG_concat.csv", "n15000_topics_plants_euNOG_concat.csv", "n20000_topics_plants_euNOG_concat.csv")
ggplot(data=df, aes(metric, y=value, group=condition, fill=condition)) +
geom_bar(stat="identity", position="dodge") +
theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
scale_fill_manual(values = c("#0072B2","#E69F00","#009E24","#FF0000", "#979797","#5530AA", "#000000"))


