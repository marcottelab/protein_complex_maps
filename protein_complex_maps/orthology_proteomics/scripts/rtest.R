
#require(ggplot2)
#require(cowplot)
#require(gdata)
#require(gplots)

print("reading data")
df <- read.csv("arathtraesorysj_euNOG_concat.csv", sep=",", header=TRUE, row.names=1)

print("clustering")
clus1 <- hclust(dist(df))

print("formatting")
d2 <- df[clus1$order,]

d2

write.csv(d2, file = "arathtraesorsj_euNOg_rowclustered.csv")

#m <- as.matrix(df, rownames.force=TRUE)
#class(m) <- "numeric"
#m[m==""] <- 0
##m[is.na(m)] <- 0
#nr <- nrow(m)
#print("rendering heatmap")
#a <- heatmap.2(head(m), scale="none", col=c("white", "black"), cexRow=0.2/log10(nr), colv=NULL, trace="none", sepcolor="grey", sepwidth=0.01, key=FALSE, xlab="DATABASES", ylab="GENES", margins=c(15,10))
#print(a)
