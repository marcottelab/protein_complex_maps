library(Hotelling)
library(huge)

#wan1.raw <- read.table('/home/kdrew/data/protein_complex_maps/flat_matrix/Hs_wan_shuhplc94_1206.txtrot_count_uniqpeps2_FDR0010_peptide_normalized_ids_mapped.txt')
#wan1.raw <- read.table('~/data/protein_complex_maps/flat_matrix/HsCeDmMmSp_all.prot_count_uniqpeps2_FDR0010_ids_mapped_proteasome.txt')
#data1.raw <- read.table('~/data/protein_complex_maps/flat_matrix/Hs_all.prot_count_uniqpeps2_FDR0010.txt')
data1.raw <- read.table('~/data/protein_complex_maps/flat_matrix/Hs_all.prot_count_uniqpeps2_FDR0010_pdb_list2.txt')
cat("read data\n")

data1.num <- sapply(data1.raw,as.numeric)
cat("data numeric\n")
rm(data1.raw)
gc()

data1.pseudo <- data1.num + 1.0
cat("added pseudo count\n")
rm(data1.num)
gc()

data1.trans <- t(data1.pseudo)
cat("transposed matrix\n")
rm(data1.pseudo)
gc()

data1.clr <- clr(data1.trans)
cat("clr transformed \n")
rm(data1.trans)
gc()

data1.npn <- huge.npn(data1.clr)
cat("npn transformed \n")
rm(data1.clr)
gc()

data1.out <- huge(data1.npn, method="mb", nlambda=30)
cat("finished huge\n")
#rm(data1.npn)
#gc()

data1.stars <- huge.select(data1.out, criterion = "stars", stars.thresh=0.05)
cat("finished stars select\n")

resGraph <- data1.stars$refit

