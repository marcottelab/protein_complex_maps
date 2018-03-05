# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
df = pd.read_table("/stor/work/Marcotte/MS/processed/Treiber_microRNA_pulldowns/group_files/T_Treiber_240212.prot_count_mFDRpsm001.elut",index_col=0)
df_resetind = df.reset_index()
df_tidy = pd.melt(df_resetind,id_vars='index')
df_tidy_nonzero = df_tidy.query("value > 0")
df_annotations = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_20171218.tab")
df_geneid_map = df_annotations[['Entry','Cross-reference (GeneID)']]
df_geneid_map.columns = ['ACC','geneid']
geneids = [x.split(';')[0] if x==x else x for x in df_geneid_map.geneid.values]
df_geneid_map['geneid_clean'] = geneids
df_geneid_map = df_geneid_map.set_index('ACC')
df.head()
df.sum()
df.sum() > 100
df
df.sum() > 100
df[df.sum() > 100]
df.loc[df.sum() > 100]
df.sum() > 100
df[df.sum() > 100]
df.columns
df.columns[df.sum()>100]
df[df.columns[df.sum()>100]]
df_100expSC = df[df.columns[df.sum()>100]]
df_100expSC_normalized = df_100expSC/df_100expSC.sum()
df_100expSC_normalized_resetind = df_100expSC_normalized.reset_index()
df_100expSC_normalized_tidy = pd.melt(df_100expSC_normalized_resetind,id_vars='index')
df_100expSC_normalized_tidy_nonzero = df_100expSC_normalized_tidy.query("value > 0")
df_100expSC_normalized_tidy_nonzero.head()
df_100expSC_normalized_scaledmax = df_100expSC_normalized*df_100expSC.max().max()
df_100expSC_normalized_scaledmax_resetind = df_100expSC_normalized_scaledmax.reset_index()
df_100expSC_normalized_scaledmax_tidy = pd.melt(df_100expSC_normalized_scaledmax_resetind,id_vars='index')
df_100expSC_normalized_scaledmax_tidy_nonzero = df_100expSC_normalized_scaledmax_tidy.query("value > 0")
df_100expSC_normalized_scaledmax_tidy_nonzero.head()
df_100expSC_normalized_scaledmax_tidy_nonzero.head(10)
df_tidy_nonzero.head()
df_tidy_nonzero.head(10)
df_merge = df_tidy_nonzero.merge(df_100expSC_normalized_scaledmax_tidy_nonzero, on=['index','variable'])
df_merge.head(10)
df_merge.head(20)
df_geneid_map
df_merge = df_merge.join(df_geneid_map, on='index')
df_merge['dataset'] = 'Treiber'
df_merge.head()
pzscores = (stats.rankdata(df_merge['value_x'])-1) / (len(df_merge['value_x'])-1) * 100
df_merge['percentile'] = pzscores
df_merge
df_merge.head()
df_merge_trim = df_merge[['variable','geneid_clean','value_x','value_y','percentile','dataset']]
df_merge_trim.head()
df_merge_trim.columns = ['experiment_id','geneid','psm','normalized_psm','percentile','dataset']
df_merge_trim.to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_nonzero_normalized.featmat")
train_df = pd.read_table("/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/complex_map2.1/corum/coreComplexes_20170702_geneids_human_merged06_rmLrg_rmSub.train_ppis.txt",header=None)
train_ppis = set(train_df[0].values).union(set(train_df[1].values))
train_ppis_str = set([str(int(x)) for x in train_ppis])
df_merge_trim.query("geneid in @train_ppis_str").to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_nonzero_normalized_train_ppis.featmat")
