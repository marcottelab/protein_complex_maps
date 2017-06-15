library(tidyr)
library(dplyr)


#This script turns a wide elution profile into a three column tidy format
#It filter to spectral counts over a certain count and removes the speccount column 
#Setup for shared_bait_feature.py

args = commandArgs(trailingOnly=TRUE)

wide_form <- read.csv(args[1], header=TRUE)

long_form <- gather(wide_form, Fraction, SpectralCounts, 2:ncol(wide_form))


#make this a function

long_form1 <- long_form %>% filter(SpectralCounts > 0 )
long_form2 <- long_form %>% filter(SpectralCounts > 1 )
long_form3 <- long_form %>% filter(SpectralCounts > 2 )
long_form4 <- long_form %>% filter(SpectralCounts > 3 )
long_form5 <- long_form %>% filter(SpectralCounts > 4 )
long_form6 <- long_form %>% filter(SpectralCounts > 5 )
long_form7 <- long_form %>% filter(SpectralCounts > 6 )
long_form8 <- long_form %>% filter(SpectralCounts > 7 )
long_form9 <- long_form %>% filter(SpectralCounts > 8 )
long_form10 <- long_form %>% filter(SpectralCounts > 9 )




long_form1 <- long_form1 %>% select(ID, Fraction)
long_form2 <- long_form2 %>% select(ID, Fraction)
long_form3 <- long_form3 %>% select(ID, Fraction)
long_form4 <- long_form4 %>% select(ID, Fraction)
long_form5 <- long_form5 %>% select(ID, Fraction)
long_form6 <- long_form6 %>% select(ID, Fraction)
long_form7 <- long_form7 %>% select(ID, Fraction)
long_form8 <- long_form8 %>% select(ID, Fraction)
long_form9 <- long_form9 %>% select(ID, Fraction)
long_form10 <- long_form10 %>% select(ID, Fraction)



write.table(long_form1, paste("min1",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form2, paste("min2",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form3, paste("min3",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form4, paste("min4",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form5, paste("min5",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form6, paste("min6",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form7, paste("min7",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form8, paste("min8",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form9, paste("min9",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")
write.table(long_form10, paste("min10",args[2], sep="_"), row.names=FALSE, quote=FALSE, sep=",")


