# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
bioplex2_fname = "/stor/work/Marcotte/project/kdrew/data/bioplex/bioplex2/BaitPreyPairs_noFilters_BP2a.tsv"
bioplex2_df = pd.read_csv(bioplex2_fname, sep='\t')
#kdrew: pandas is a bit annoying sometimes when dealing with numbers that should be strings but are considered floats
#kdrew: get uniprot ids of entries where gene_id is NaN
nan_prot_id = list(bioplex2_df.query("symbol != symbol")['db_protein_id'].values)
#kdrew: set gene_id NaNs to be an sentinal int (-1)
bioplex2_df.loc[bioplex2_df["symbol"].isnull(),'gene_id'] = -1
#kdrew: convert all gene_ids into ints and then strings
bioplex2_df['gene_id_str'] = bioplex2_df['gene_id'].apply(int).apply(str)
#kdrew: replace sentinal with uniprot ids for entries that were originally NaN
bioplex2_df.loc[bioplex2_df["symbol"].isnull(),'gene_id_str'] = nan_prot_id
orig9k_fname = "/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_feature_revisitMerge_pairsOnly_preyMerge2_heinCollapseMerge_pairsOnly_preyMerge2.txt"
orig9k_df = pd.read_csv(orig9k_fname)
get_ipython().magic(u'save humap2_create_featmat 1-20')
orig9k_df.head()
orig9k_df.columns
orig9k_df.columns.tolist()
orig9k_df['frozenset_geneids']
orig9k_df['frozenset_geneids_str_order']
orig9k_df['frozenset_geneids_str_order'].isnan()
orig9k_df.query("frozenset_geneids_str_order != frozenset_geneids_str_order")
orig9k_df.query("frozenset_geneids_str_order == frozenset_geneids_str_order")
orig9k_df['frozenset_geneids_str_order']
orig9k_df['frozenset_geneids_str_order'].values[0]
bioplex2_df.head()
orig9k_df.columns
orig9k_df.columns.tolist()
orig9k_df.columns.tolist()[3:]
orig9k_df.columns.tolist()[3:200]
orig9k_df.columns.tolist()[3:205]
orig9k_df.columns.tolist()[3:210]
orig9k_df.columns.tolist()[3:220]
orig9k_df.columns.tolist()[3:225]
orig9k_df.columns.tolist()[3:223]
orig9k_df.columns.tolist()[3:223,230]
orig9k_df.columns.tolist()[3:223]
wan_frac_feats = orig9k_df.columns.tolist()[3:223]
orig9k_df.columns.tolist()[223:]
orig9k_df.columns.tolist()[222:]
orig9k_df.columns.tolist()[230:]
orig9k_df.columns.tolist()[240:]
orig9k_df.columns.tolist()[242:]
orig9k_df.columns.tolist()[242:244]
wan_extra_feat = orig9k_df.columns.tolist()[242:244]
orig9k_df.columns.tolist()[244:]
orig9k_df.columns.tolist()[246:]
orig9k_df.columns.tolist()[247:]
orig9k_df.columns.tolist()[249:]
orig9k_df.columns.tolist()[249:257]
orig9k_df.columns.tolist()[249:258]
bioplex_feats = orig9k_df.columns.tolist()[249:258]
bioplex_feats = orig9k_df.columns.tolist()[258:]
bioplex_feats = orig9k_df.columns.tolist()[249:258]
orig9k__df.columns.tolist()[258:]
orig9k_df.columns.tolist()[258:]
orig9k_df.columns.tolist()[260:]
orig9k_df.columns.tolist()[260:262]
bioplex_prey_prey_feats = orig9k_df.columns.tolist()[260:262]
orig9k_df.columns.tolist()[262:]
orig9k_df.columns.tolist()[265:]
orig9k_df.columns.tolist()[265:270]
hein_feats = orig9k_df.columns.tolist()[265:270]
orig9k_df.columns.tolist()[270:]
orig9k_df.columns.tolist()[272:]
orig9k_df.columns.tolist()[272:274]
hein_prey_prey_feats = orig9k_df.columns.tolist()[272:274]
get_ipython().magic(u'lsmagic')
vars
vars()
feat_names = dict()
feat_names['wan_frac_feats'] = wan_frac_feats
feat_names['wan_extra_feat'] = wan_extra_feat
feat_names['bioplex_feats'] = bioplex_feats
feat_names['bioplex_prey_prey_feats'] = bioplex_prey_prey_feats
feat_names['hein_feats'] = hein_feats
feat_names['hein_prey_prey_feats'] = hein_prey_prey_feats
get_ipython().magic(u'save humap2_create_featmat 1-83')
feat_names
orig9k_df.columns.tolist()
orig9k_trim_columns = ["frozenset_geneids_str_order"] + feat_names['wan_frac_feats'] + feat_names['wan_extra_feat'] + feat_names['bioplex_feats'] + feat_names['bioplex_prey_prey_feats'] + feat_names['hein_feats'] + feat_names['hein_prey_prey_feats']
orig9k_trim_columns
orig9k_df[orig9k_trim_columns]
orig9k_trim_columns = ['geneid1', 'geneid2'] + orig9k_trim_columns
orig9k_trim_columns
orig9k_df[orig9k_trim_columns].head()
orig9k_trim_df = orig9k_df[orig9k_trim_columns]
orig9k_trim_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k.featmat")
bioplex2_df.head()
bioplex2_df.columns.tolist()
bioplex2_df[['bait_geneid','gene_id_str','ave_apsm','nwdscore','zscore','plate_zscore','entropy','uPeps','ratio','total_psms','ratioTotalPSMs','UtoTratio']].head()
bioplex2_trim_df = bioplex2_df[['bait_geneid','gene_id_str','ave_apsm','nwdscore','zscore','plate_zscore','entropy','uPeps','ratio','total_psms','ratioTotalPSMs','UtoTratio']]
bioplex2_trim_df.to_csv("/stor/work/Marcotte/project/kdrew/data/bioplex/bioplex2/bioplex2.featmat")
bioplex2_trim_df.to_csv("/stor/work/Marcotte/project/kdrew/data/bioplex/bioplex2/bioplex2.featmat", index=False)
orig9k_trim_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k.featmat", index=False)
