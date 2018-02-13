# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
bioid_df = pd.read_csv("/project/kdrew/data/youn_bioid/youn_bioid_mmc2-5_SAINT_ID3105_wgeneids.csv")
bioid_df.head()
youn_bioid_featmat_df = bioid_df[['Bait Name','Bait Gene','Prey Gene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','geneids','prey_gene_id']]
youn_bioid_featmat_df = bioid_df[['Bait Name','Bait Gene','Prey Gene','AvgSpec','AvgP','MaxP','FoldChange','BFDR','gene_ids','prey_gene_ids']]
youn_bioid_featmat_df['geneid_pairs_strsort'] = [str(sorted( [ str(x[0]), str(x[1]) ] )) for x in youn_bioid_featmat_df[['gene_ids','prey_gene_ids']].values] 
youn_bioid_featmat_df = youn_bioid_featmat_df.set_index("geneid_pairs_strsort")
youn_bioid_featmat_df.to_csv("/stor/work/Marcotte/project/kdrew/data/youn_bioid/youn_bioid_All.featmat")
youn_bioid_hypergeo_input_df = youn_bioid_featmat_df.reset_index()[['Bait Name','prey_gene_ids','AvgSpec']]
youn_bioid_hypergeo_input_df['dataset'] = 'youn'
pzscores = (stats.rankdata(youn_bioid_hypergeo_input_df['AvgSpec'])-1) / (len(youn_bioid_hypergeo_input_df['AvgSpec'])-1) * 100
youn_bioid_hypergeo_input_df['abundance'] = pzscores
youn_bioid_hypergeo_input_df.head()
youn_hypergeo_bioidonly_df.columns = ['experiment_id','geneid','AvgSpec','dataset','abundance']
youn_hypergeo_input_df.columns = ['experiment_id','geneid','AvgSpec','dataset','abundance']
youn_bioid_hypergeo_input_df.columns = ['experiment_id','geneid','AvgSpec','dataset','abundance']
youn_bioid_hypergeo_df[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/youn_bioid/youn_bioid_hypergeo_input.featmat")
youn_bioid_hypergeo_input_df[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/youn_bioid/youn_bioid_hypergeo_input.featmat")
youn_bioid_hypergeo_input_df.query("AvgSpec >= 2")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/youn_bioid/youn_bioid_hypergeo_input_AvgSpec2.featmat")
youn_bioid_hypergeo_input_df.query("AvgSpec >= 4")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/youn_bioid/youn_bioid_hypergeo_input_AvgSpec4.featmat")
