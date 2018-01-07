# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
orig9k_trim_df = pd.read_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k.featmat")
bioplex2_trim_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/bioplex/bioplex2/bioplex2.featmat")
bioplex2_geneid_pairs = bioplex2_trim_df[['bait_geneid','gene_id_str']].values
bioplex2_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in bioplex2_geneid_pairs]
bioplex2_trim_df['geneids_str_order'] = bioplex2_geneid_pairs_strsort
bioplex2_trim_df.head()
bioplex2_trim_df.set_index("geneids_str_order")
bioplex2_trim_df = bioplex2_trim_df.set_index("geneids_str_order") 
orig9k_trim_df.head()
orig9k_trim_df.head()
orig9k_trim_df = orig9k_trim_df.set_index("frozenset_geneids_str_order")
orig9k_trim_df.head()
orig9k_trim_df.head(2)
bioplex2_trim_df.head()
orig9k_bioplex2_df = orig9k_trim_df.join(bioplex2_trim_df, how="outer")
orig9k_bioplex2_df = orig9k_trim_df.join(bioplex2_trim_df, how="outer", lsuffix="_orig9k", rsuffix="_bioplex2")
orig9k_bioplex2_df.head()
orig9k_bioplex2_df.tail()
orig9k_bioplex2_df.columns
orig9k_bioplex2_df.query("geneid1 == geneid1 and bait_geneid == bait_geneid").head()
orig9k_bioplex2_df.query("Ce_6mg_1203_pq_euc == Ce_6mg_1203_pq_euc and ratioTotalPSMs_bioplex2 == ratioTotalPSMs_bioplex2").head()
orig9k_bioplex2_df.query("Ce_6mg_1203_pq_euc == Ce_6mg_1203_pq_euc and ratioTotalPSMs_bioplex2 == ratioTotalPSMs_bioplex2").count()
get_ipython().magic(u'save humap2_create_featmat2 1:30')
orig9k_bioplex2_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2.featmat")
bioplex2_prey_Z4_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/combined_dataframes/BaitPreyPairs_noFilters_BP2a_Z4_prey_pairs_pvalcorr_logchoose")
bioplex2_prey_Z2_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/combined_dataframes/BaitPreyPairs_noFilters_BP2a_Z2_prey_pairs_pvalcorr_logchoose")
bioplex2_prey_Z4_df.head()
bioplex2_prey_Z4_geneid_pairs = bioplex2_prey_Z4_df[['gene_id1','gene_id2']].values
bioplex2_prey_Z2_geneid_pairs = bioplex2_prey_Z2_df[['gene_id1','gene_id2']].values
bioplex2_prey_Z4_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in bioplex2_prey_Z4_geneid_pairs]
bioplex2_prey_Z2_geneid_pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in bioplex2_prey_Z2_geneid_pairs]
bioplex2_prey_Z4_df['geneids_str_order'] = bioplex2_prey_Z4_geneid_pairs_strsort
bioplex2_prey_Z2_df['geneids_str_order'] = bioplex2_prey_Z2_geneid_pairs_strsort
bioplex2_prey_Z4_df = bioplex2_prey_Z4_df.set_index("geneids_str_order")
bioplex2_prey_Z2_df = bioplex2_prey_Z2_df.set_index("geneids_str_order")
orig9k_bioplex2_hygeoZ4 = orig9k_bioplex2_df.join(bioplex2_prey_Z4_df, how="outer", lsuffix="_orig9k", rsuffix="_bioplex2_Z4")
orig9k_bioplex2_hygeoZ4_Z2_df = orig9k_bioplex2_hygeoZ4.join(bioplex2_prey_Z2_df, how="outer", lsuffix="_orig9k", rsuffix="_bioplex2_Z2")
orig9k_bioplex2_hygeoZ4_Z2_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2.featmat")
get_ipython().magic(u'save humap2_create_featmat2 1:46')
orig9k_bioplex2_hygeoZ4 = orig9k_bioplex2_df.join(bioplex2_prey_Z4_df, how="outer", rsuffix="_bioplex2_Z4")
orig9k_bioplex2_hygeoZ4_Z2_df = orig9k_bioplex2_hygeoZ4.join(bioplex2_prey_Z2_df, how="outer", rsuffix="_bioplex2_Z2")
orig9k_bioplex2_hygeoZ4_Z2_df.columns.tolist()
orig9k_bioplex2_hygeoZ4_Z2_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2.featmat")
orig9k_bioplex2_hygeoZ4_Z2_df.head()
orig9k_bioplex2_hygeoZ4_Z2_df.query("Ce6mg_1203_pq_euc == Ce6mg_1203_pq_euc and pval_bioplex2_Z2 == pval_bioplex2_Z2").head()
orig9k_bioplex2_hygeoZ4_Z2_df.head()
orig9k_bioplex2_hygeoZ4_Z2_df.query("Ce6mg_1203_apex == Ce6mg_1203_apex and pval_bioplex2_Z2 == pval_bioplex2_Z2").head()
orig9k_bioplex2_hygeoZ4_Z2_df.head()
orig9k_bioplex2_hygeoZ4_Z2_df.columns()
orig9k_bioplex2_hygeoZ4_Z2_df.columns
orig9k_bioplex2_hygeoZ4_Z2_df.query("Ce_6mg_1203_apex == Ce_6mg_1203_apex and pval_bioplex2_Z2 == pval_bioplex2_Z2").head()
orig9k_bioplex2_hygeoZ4_Z2_df.query("Ce_6mg_1203_apex == Ce_6mg_1203_apex and pval_bioplex2_Z2 == pval_bioplex2_Z2").count()
get_ipython().magic(u'save humap2_create_featmat2 1:60')
orig9k_bioplex2_hygeoZ4_Z2_df.head()
ciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_ciliated.featmat")
ciliated_bioid_df.head()
nonciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_non_ciliated.featmat")
nonciliated_bioid_df.head()
nonciliated_bioid_df = nonciliated_bioid_df.set_index("geneid_pairs_strsort")
ciliated_bioid_df = ciliated_bioid_df.set_index("geneid_pairs_strsort")
ciliated_bioid_df.head()
orig9k_bioplex2_hygeoZ4_Z2_cbioid_df = orig9k_bioplex2_hygeoZ4_Z2_df.join(ciliated_bioid_df, how="outer", rsuffix="_ciliated_bioid")
orig9k_bioplex2_hygeoZ4_Z2_cbioid_df.head()
orig9k_bioplex2_hygeoZ4_Z2_cbioid_df.query("Ce_6mg_1203_poisson == Ce_6mg_1203_poisson and AvgP == AvgP").head()
orig9k_bioplex2_hygeoZ4_Z2_bioid_df = orig9k_bioplex2_hygeoZ4_Z2_cbioid_df.join(nonciliated_bioid_df, how="outer", rsuffix="_nonciliated_bioid")
orig9k_bioplex2_hygeoZ4_Z2_bioid_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2_bioid.featmat")
bioplex2_prey_Z2_gb_df = bioplex2_prey_Z2_df.groupby(bioplex2_prey_Z2_df.index).first()
bioplex2_prey_Z4_gb_df = bioplex2_prey_Z4_df.groupby(bioplex2_prey_Z4_df.index).first()
len(bioplex2_prey_Z2_gb_df)
len(bioplex2_prey_Z2_df)
len(bioplex2_prey_Z4_gb_df)
len(bioplex2_prey_Z4_df)
orig9k_bioplex2_hygeoZ4_df = orig9k_bioplex2_df.join(bioplex2_prey_Z4_gb_df, how="outer", rsuffix="_bioplex2_Z4")
orig9k_bioplex2_hygeoZ4_Z2_df = orig9k_bioplex2_hygeoZ4_df.join(bioplex2_prey_Z2_gb_df, how="outer", rsuffix="_bioplex2_Z2")
len(orig9k_bioplex2_hygeoZ4_Z2_df)
orig9k_bioplex2_hygeoZ4_Z2_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2.featmat")
