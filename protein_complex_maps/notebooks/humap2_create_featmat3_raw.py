# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
orig9k_bioplex2_hygeo_df = pd.read_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2.featmat")
orig9k_bioplex2_hygeo_df.head()
orig9k_bioplex2_hygeo_df.head().values
orig9k_bioplex2_hygeo_df.head().values[0]
orig9k_bioplex2_hygeo_df.head().values[1]
orig9k_bioplex2_hygeo_df.head().values[10]
orig9k_bioplex2_hygeo_df.head().values[3]
orig9k_bioplex2_hygeo_df.values[30]
orig9k_bioplex2_hygeo_df.ix[30]
orig9k_bioplex2_hygeo_df.ix[3]
orig9k_bioplex2_hygeo_df.ix[300]
orig9k_bioplex2_hygeo_df.ix[3000]
orig9k_bioplex2_hygeo_df.ix[30000]
orig9k_bioplex2_hygeo_df.ix[300000]
orig9k_bioplex2_hygeo_df.ix[3000000]
orig9k_bioplex2_hygeo_df.ix[30000000]
orig9k_bioplex2_hygeo_df.ix[3000000]
orig9k_bioplex2_hygeo_df.ix[900]
orig9k_bioplex2_hygeo_df.ix[9000]
orig9k_bioplex2_hygeo_df.ix[90000]
orig9k_bioplex2_hygeo_df.ix[900000]
orig9k_bioplex2_hygeo_df.ix[90000]
orig9k_bioplex2_hygeo_df.ix[80000]
orig9k_bioplex2_hygeo_df.ix[70000]
orig9k_bioplex2_hygeo_df.ix[40000]
ciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_ciliated.featmat")
nonciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_non_ciliated.featmat")
nonciliated_bioid_df = nonciliated_bioid_df.set_index("geneid_pairs_strsort")
ciliated_bioid_df = ciliated_bioid_df.set_index("geneid_pairs_strsort")
ciliated_bioid_df.head()
len(ciliated_bioid_df)
len(ciliated_bioid_df.groupby(ciliated_bioid_df.index).first())
ciliated_bioid_gb_df = ciliated_bioid_df.groupby(ciliated_bioid_df.index).first()
nonciliated_bioid_gb_df = nonciliated_bioid_df.groupby(nonciliated_bioid_df.index).first()
len(nonciliated_bioid_gb_df)
len(nonciliated_bioid_df)
ciliated_bioid_df.groupby(ciliated_bioid_df.index).agg(['count'])
ciliated_bioid_df.index.groupby().agg(['count'])
ciliated_bioid_df.index.groupby(ciliated_bioid_df.index).agg(['count'])
ciliated_bioid_df.index.groupby(ciliated_bioid_df.index)
ciliated_bioid_df.groupby(ciliated_bioid_df.index).size()
ciliated_bioid_df.groupby(ciliated_bioid_df.index).size().sort()
ciliated_bioid_df.groupby(ciliated_bioid_df.index).size().sort_values()
ciliated_bioid_df.query("index == '[\'2778\', \'80776\']'")
ciliated_bioid_df.iloc["['2778', '80776']"]
ciliated_bioid_df[ciliated_bioid_df.index == "['2778', '80776']"]
ciliated_bioid_df[ciliated_bioid_df.index == "['11064', '2778']"]
ciliated_bioid_df[ciliated_bioid_df.index == "['7112', '80184']"]
orig9k_bioplex2_hygeo_df.head()
orig9k_bioplex2_hygeo_cbioid_df = orig9k_bioplex2_hygeo_df.join(ciliated_bioid_gb_df, how="outer", rsuffix="_ciliated_bioid")
orig9k_bioplex2_hygeo_cbioid_df.head()
orig9k_bioplex2_hygeo_df = orig9k_bioplex2_hygeo_df.set_index("Unnamed: 0")
orig9k_bioplex2_hygeo_df.head()
orig9k_bioplex2_hygeo_cbioid_df = orig9k_bioplex2_hygeo_df.join(ciliated_bioid_gb_df, how="outer", rsuffix="_ciliated_bioid")
orig9k_bioplex2_hygeo_bioid_df = orig9k_bioplex2_hygeo_cbioid_df.join(nonciliated_bioid_gb_df, how="outer", rsuffix="_nonciliated_bioid")
orig9k_bioplex2_hygeo_bioid_df.head()
orig9k_bioplex2_hygeo_bioid_df.query("Ce_6mg_1203_pq_euc == Ce_6mg_1203_pq_euc and condition_nonciliated_bioid == condition_nonciliated_bioid").head()
orig9k_bioplex2_hygeo_bioid_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeo_bioid.featmat")
len(orig9k_bioplex2_hygeo_bioid_df)
orig9k_bioplex2_hygeo_bioid_df.head()
get_ipython().magic(u'save humap2_create_featmat3 1:68')
cilium_hygeo_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_prey_pairs_pvalcorr_logchoose.feat")
cilium_hygeo_df.head()
cilium_hygeo_avgspec2_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_AvgSpec2_prey_pairs_pvalcorr_logchoose.feat") 
cilium_hygeo_avgspec4_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_AvgSpec4_prey_pairs_pvalcorr_logchoose.feat")
cilium_hygeo_geneid_pairs = cilium_hygeo_df[['gene_id1','gene_id2']].values
cilium_hygeo_avgspec2_geneid_pairs = cilium_hygeo_avgspec2_df[['gene_id1','gene_id2']].values
cilium_hygeo_avgspec4_geneid_pairs = cilium_hygeo_avgspec4_df[['gene_id1','gene_id2']].values
cilium_hygeo_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in cilium_hygeo_geneid_pairs]
cilium_hygeo_avgspec2_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in cilium_hygeo_avgspec2_geneid_pairs]
cilium_hygeo_avgspec4_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in cilium_hygeo_avgspec4_geneid_pairs]
cilium_hygeo_df['geneids_str_order'] = cilium_hygeo_geneid_pairs_strsort
cilium_hygeo_avgspec2_df['geneids_str_order'] = cilium_hygeo_avgspec2_geneid_pairs_strsort
cilium_hygeo_avgspec4_df['geneids_str_order'] = cilium_hygeo_avgspec4_geneid_pairs_strsort
cilium_hygeo_df = cilium_hygeo_df.set_index("geneids_str_order")
cilium_hygeo_avgspec2_df = cilium_hygeo_avgspec2_df.set_index("geneids_str_order")
cilium_hygeo_avgspec4_df = cilium_hygeo_avgspec4_df.set_index("geneids_str_order")
bioplex2_prey_Z2_gb_df = bioplex2_prey_Z2_df.groupby(bioplex2_prey_Z2_df.index).first()
cilium_hygeo_gb_df = cilium_hygeo_df.groupby(cilium_hygeo_df.index).first()
cilium_hygeo_avgspec2_gb_df = cilium_hygeo_avgspec2_df.groupby(cilium_hygeo_avgspec2_df.index).first()
cilium_hygeo_avgspec4_gb_df = cilium_hygeo_avgspec4_df.groupby(cilium_hygeo_avgspec4_df.index).first()
orig9k_bioplex2_hygeo_bioid_hygeo_df = orig9k_bioplex2_hygeo_bioid_df.join(cilium_hygeo_gb_df, how="outer", rsuffix="_cilium_hygeo")
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_df = orig9k_bioplex2_hygeo_bioid_hygeo_df.join(cilium_hygeo_avgspec2_gb_df, how="outer", rsuffix="_cilium_hygeo_avgspec2")
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df = orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_df.join(cilium_hygeo_avgspec4_gb_df, how="outer", rsuffix="_cilium_hygeo_avgspec4")
len(orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df)
20000*20000
1.0*len(orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df) /400000000
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df.head()
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeo_bioid_hygeo.featmat")
