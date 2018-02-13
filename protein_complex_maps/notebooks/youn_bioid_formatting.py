# coding: utf-8
#kdrew: all code was run on pocketfox
import pandas as pd

df2 = pd.read_csv("youn_bioid_mmc2-5_SAINT_ID3105.csv")

f = open("./prey_names_20180213.txt","wb")
f.write("\n".join(df2['Prey Gene'].unique()))
f.close()

#kdrew: search uniprot with genenames
#kdrew: searched bait genenames on uniprot manually
bait_uniprot_df = pd.read_table("./uniprot-bait_names_20180213.tab")

bait_geneids = [x.split(';')[0] for x in bait_uniprot_df['Cross-reference (GeneID)'].values]
bait_uniprot_df['gene_ids'] = bait_geneids

bait_uniprot_df = bait_uniprot_df.set_index('yourlist:M20180213AAFB7E4D2F1D05654627429E83DA5CCEED911CP')
df2_bait_uniprot = df2.join(bait_uniprot_df['gene_ids'],on='Bait Gene')

#kdrew: search uniprot with genenames
prey_uniprot_df = pd.read_table("./uniprot-prey_names_20180213.tab")

prey_geneids = [str(x).split(';')[0] for x in prey_uniprot_df['Cross-reference (GeneID)'].values]
prey_uniprot_df['prey_gene_ids'] = prey_geneids

prey_uniprot_df = prey_uniprot_df.set_index('yourlist:M20180213A7434721E10EE6586998A056CCD0537E3D64557')
df2_bait_prey_uniprot = df2_bait_uniprot.join(prey_uniprot_df['prey_gene_ids'], on='Prey Gene')

df2_bait_prey_uniprot.to_csv("./youn_bioid_mmc2-5_SAINT_ID3105_wgeneids.csv",index=None)

