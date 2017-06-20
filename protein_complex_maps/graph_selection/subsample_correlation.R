
library(Hotelling)
library(huge)

data1.raw <- read.table('/project/kdrew/data/protein_complex_maps/flat_matrix/Hs_all.prot_count_uniqpeps2_FDR0010.txt')
cat("read data\n")

data1.num <- sapply(data1.raw,as.numeric)
cat("data numeric\n")
rm(data1.raw)
gc()


columns <- colnames(data1.num)
shuffled_columns <- sample(columns) 

d1 <- data1.num[,shuffled_columns[1:150]]
d2 <- data1.num[,shuffled_columns[151:300]]
d3 <- data1.num[,shuffled_columns[301:450]]
d4 <- data1.num[,shuffled_columns[451:600]]
d5 <- data1.num[,shuffled_columns[601:750]]
d6 <- data1.num[,shuffled_columns[751:900]]
d7 <- data1.num[,shuffled_columns[901:1050]]
d8 <- data1.num[,shuffled_columns[1051:1200]]
d9 <- data1.num[,shuffled_columns[1201:1350]]
d10 <- data1.num[,shuffled_columns[1351:1500]]
d11 <- data1.num[,shuffled_columns[1501:1650]]
d12 <- data1.num[,shuffled_columns[1651:1800]]
d13 <- data1.num[,shuffled_columns[1801:1950]]
d14 <- data1.num[,shuffled_columns[1951:2100]]
d15 <- data1.num[,shuffled_columns[2101:2250]]
d16 <- data1.num[,shuffled_columns[2251:2400]]
d17 <- data1.num[,shuffled_columns[2401:2550]]
d18 <- data1.num[,shuffled_columns[2551:2700]]
d19 <- data1.num[,shuffled_columns[2701:2850]]
d20 <- data1.num[,shuffled_columns[2851:2989]]

cor_d1 <- cor(t(d1))
cor_d1_0 <- replace(cor_d1, is.na(cor_d1), 0)
rm(d1)
rm(cor_d1)
gc()

cor_d2 <- cor(t(d2))
cor_d2_0 <- replace(cor_d2, is.na(cor_d2), 0)
rm(d2)
rm(cor_d2)
gc()

cor_d3 <- cor(t(d3))
cor_d3_0 <- replace(cor_d3, is.na(cor_d3), 0)
rm(d3)
rm(cor_d3)
gc()

cor_d4 <- cor(t(d4))
cor_d4_0 <- replace(cor_d4, is.na(cor_d4), 0)
rm(d4)
rm(cor_d4)
gc()

cor_d5 <- cor(t(d5))
cor_d5_0 <- replace(cor_d5, is.na(cor_d5), 0)
rm(d5)
rm(cor_d5)
gc()

cor_d6 <- cor(t(d6))
cor_d6_0 <- replace(cor_d6, is.na(cor_d6), 0)
rm(d6)
rm(cor_d6)
gc()

cor_d7 <- cor(t(d7))
cor_d7_0 <- replace(cor_d7, is.na(cor_d7), 0)
rm(d7)
rm(cor_d7)
gc()

cor_d8 <- cor(t(d8))
cor_d8_0 <- replace(cor_d8, is.na(cor_d8), 0)
rm(d8)
rm(cor_d8)
gc()

cor_d9 <- cor(t(d9))
cor_d9_0 <- replace(cor_d9, is.na(cor_d9), 0)
rm(d9)
rm(cor_d9)
gc()

cor_d10 <- cor(t(d10))
cor_d10_0 <- replace(cor_d10, is.na(cor_d10), 0)
rm(d10)
rm(cor_d10)
gc()

cor_d11 <- cor(t(d11))
cor_d11_0 <- replace(cor_d11, is.na(cor_d11), 0)
rm(d11)
rm(cor_d11)
gc()

cor_d12 <- cor(t(d12))
cor_d12_0 <- replace(cor_d12, is.na(cor_d12), 0)
rm(d12)
rm(cor_d12)
gc()

cor_d13 <- cor(t(d13))
cor_d13_0 <- replace(cor_d13, is.na(cor_d13), 0)
rm(d13)
rm(cor_d13)
gc()

cor_d14 <- cor(t(d14))
cor_d14_0 <- replace(cor_d14, is.na(cor_d14), 0)
rm(d14)
rm(cor_d14)
gc()

cor_d15 <- cor(t(d15))
cor_d15_0 <- replace(cor_d15, is.na(cor_d15), 0)
rm(d15)
rm(cor_d15)
gc()

cor_d16 <- cor(t(d16))
cor_d16_0 <- replace(cor_d16, is.na(cor_d16), 0)
rm(d16)
rm(cor_d16)
gc()

cor_d17 <- cor(t(d17))
cor_d17_0 <- replace(cor_d17, is.na(cor_d17), 0)
rm(d17)
rm(cor_d17)
gc()

cor_d18 <- cor(t(d18))
cor_d18_0 <- replace(cor_d18, is.na(cor_d18), 0)
rm(d18)
rm(cor_d18)
gc()

cor_d19 <- cor(t(d19))
cor_d19_0 <- replace(cor_d19, is.na(cor_d19), 0)
rm(d19)
rm(cor_d19)
gc()

cor_d20 <- cor(t(d20))
cor_d20_0 <- replace(cor_d20, is.na(cor_d20), 0)
rm(d20)
rm(cor_d20)
gc()

cor_dsum <- cor_d1_0 + cor_d2_0 + cor_d3_0 + cor_d4_0 + cor_d5_0 + cor_d6_0 + cor_d7_0 + cor_d8_0 + cor_d9_0 + cor_d10_0 + cor_d11_0 + cor_d12_0 + cor_d13_0 + cor_d14_0 + cor_d15_0 + cor_d16_0 + cor_d17_0 + cor_d18_0 + cor_d19_0 + cor_d20_0

cor_davg <- cor_dsum/20

write.table(cor_davg, file="/project/kdrew/data/protein_complex_maps/graph_selection/Hs_all/correlations/Hs_all.prot_count_uniqpeps2_FDR0010.correlation.subsampled20.txt", row.names=FALSE, col.names=FALSE)



