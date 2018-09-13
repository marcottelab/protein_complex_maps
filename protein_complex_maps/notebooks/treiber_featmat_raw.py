# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
df = pd.read_table("/stor/work/Marcotte/MS/processed/Treiber_microRNA_pulldowns/group_files/T_Treiber_240212.prot_count_mFDRpsm001.elut")
df.head()
df = pd.read_table("/stor/work/Marcotte/MS/processed/Treiber_microRNA_pulldowns/group_files/T_Treiber_240212.prot_count_mFDRpsm001.elut",index_col=0)
df.head()
df_tidy = df.melt()
df_tidy = pd.melt(df)
df_tidy.head()
df_tidy.head(100)
df_melt.head()
df_tidy.head()
len(df_tidy)
df_nonzero_tidy = df_tidy.query("value > 0")
len(df_nonzero_tidy )
df_nonzero_tidy.head()
df_nonzero_tidy.head(100)
df.head()
df.reset_index()
df_resetind = df.reset_index()
df_resetind.head()
df_tidy = pd.melt(df,id_vars='index')
df_tidy = pd.melt(df_resetind,id_vars='index')
df_tidy.head()
df_tidy_nonzero = df_tidy.query("value > 0")
df_tidy_nonzero.head()
df_tidy_nonzero.head(100)
df_tidy_nonzero.head()
from scipy import stats
pzscores = (stats.rankdata(df_tidy_nonzero['value'])-1) / (len(df_tidy_nonzero['value'])-1) * 100
df_tidy_nonzero['abundance'] = pzscores
df_tidy_nonzero.head()
df_tidy_gt2 = df_tidy.query("value > 2")
pzscores_gt2 = (stats.rankdata(df_tidy_gt2['value'])-1) / (len(df_tidy_gt2['value'])-1) * 100
df_tidy_gt2['abundance'] = pzscores_gt2
df_tidy_gt2.head()
df_annotations = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_20171218.tab")
df_annotations.head()
df_annotations[['Entry','Cross-reference (GeneID)']].head()
df_geneid_map = df_annotations[['Entry','Cross-reference (GeneID)']]
df_geneid_map.columns = ['ACC','geneid']
df_geneid_map.head()
df_geneid_map.geneid.head()
df_geneid_map.geneid.values
geneids = [x.split(';')[0] for x in df_geneid_map.geneid.values if not np.isnan(x) else x]
geneids = [x.split(';')[0] if not np.isnan(x) else x for x in df_geneid_map.geneid.values]
df_geneid_map.geneid.head()
df_geneid_map.geneid.values[4]
np.isnan(df_geneid_map.geneid.values[4])
np.isnan(df_geneid_map.geneid.values[2])
df_geneid_map.geneid.values[4]
df_geneid_map.geneid.values[4] = df_geneid_map.geneid.values[4]
df_geneid_map.geneid.values[4] == df_geneid_map.geneid.values[4]
df_geneid_map.geneid.values[2] == df_geneid_map.geneid.values[2]
geneids = [x.split(';')[0] if x==x else x for x in df_geneid_map.geneid.values]
geneids[:10]
df_geneid_map['geneid_clean'] = geneids
df_geneid_map.head()
df_geneid_map.index('ACC')
df_geneid_map.set_index('ACC')
df_geneid_map = df_geneid_map.set_index('ACC')
df_geneid_map.head()
df_tidy_gt2.head()
df_tidy_gt2.join(df_geneid_map, index='index')
df_tidy_gt2.join(df_geneid_map, on='index')
df_tidy_gt2.join(df_geneid_map, on='index').head()
df_tidy_gt2 = df_tidy_gt2.join(df_geneid_map, on='index')
df_tidy_nonzero = df_tidy_nonzero.join(df_geneid_map, on='index')
len(df_tidy_gt2)
len(df_tidy_nonzero)
df_tidy_nonzero['dataset'] = 'Treiber'
df_tidy_gt2['dataset'] = 'Treiber'
df_tidy_gt2.head()
df_tidy_gt2_trim = df_tidy_gt2[['variable','geneid_clean','abundance','dataset']]
df_tidy_nonzero_trim = df_tidy_gt2[['variable','geneid_clean','abundance','dataset']] 
df_tidy_gt2_trim.head()
df_tidy_gt2_trim.columns = ['experiment_id','geneid','abundance','dataset']
df_tidy_nonzero_trim.columns = ['experiment_id','geneid','abundance','dataset']
df_tidy_gt2_trim.head()
df_tidy_gt2_trim = df_tidy_gt2_trim.reindex() 
df_tidy_nonzero_trim = df_tidy_nonzero_trim.reindex()
df_tidy_nonzero_trim.head()
df_tidy_nonzero_trim.to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_nonzero.featmat")
df_tidy_gt2_trim.to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_gt2.featmat")