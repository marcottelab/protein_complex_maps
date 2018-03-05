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

df_tidy_nonzero = df_tidy_nonzero.join(df_geneid_map, on='index')
df_tidy_nonzero['dataset'] = 'Treiber'
pzscores = (stats.rankdata(df_tidy_nonzero['value'])-1) / (len(df_tidy_nonzero['value'])-1) * 100
df_tidy_nonzero['percentile'] = pzscores
df_tidy_nonzero_trim = df_tidy_nonzero[['variable','geneid_clean','value','percentile','dataset']]
df_tidy_nonzero_trim.columns = ['experiment_id','geneid','psm','percentile','dataset']
df_tidy_nonzero_trim.to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_nonzero.featmat")

#get_ipython().magic(u'save treiber_featmat_20180305_raw 1:30')

#kdrew: feature matrix with only training ppis ids for testing purposes
train_df = pd.read_table("/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/complex_map2.1/corum/coreComplexes_20170702_geneids_human_merged06_rmLrg_rmSub.train_ppis.txt",header=None)
train_ppis = set(train_df[0].values).union(set(train_df[1].values))
df_tidy_nonzero_trim.query("geneid in @train_ppis").head()
train_ppis_str = set([str(int(x)) for x in train_ppis])
df_tidy_nonzero_trim.query("geneid in @train_ppis_str").to_csv("/project/kdrew/data/Treiber_etal/t_treiber_240212_nonzero_train_ppis.featmat")

