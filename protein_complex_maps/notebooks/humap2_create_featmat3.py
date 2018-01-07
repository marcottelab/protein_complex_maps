# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats

orig9k_bioplex2_hygeo_df = pd.read_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeoZ4_Z2.featmat")
orig9k_bioplex2_hygeo_df = orig9k_bioplex2_hygeo_df.set_index("Unnamed: 0")

ciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_ciliated.featmat")
nonciliated_bioid_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_non_ciliated.featmat")

nonciliated_bioid_df = nonciliated_bioid_df.set_index("geneid_pairs_strsort")
ciliated_bioid_df = ciliated_bioid_df.set_index("geneid_pairs_strsort")

ciliated_bioid_gb_df = ciliated_bioid_df.groupby(ciliated_bioid_df.index).first()
nonciliated_bioid_gb_df = nonciliated_bioid_df.groupby(nonciliated_bioid_df.index).first()

orig9k_bioplex2_hygeo_cbioid_df = orig9k_bioplex2_hygeo_df.join(ciliated_bioid_gb_df, how="outer", rsuffix="_ciliated_bioid")
orig9k_bioplex2_hygeo_bioid_df = orig9k_bioplex2_hygeo_cbioid_df.join(nonciliated_bioid_gb_df, how="outer", rsuffix="_nonciliated_bioid")

orig9k_bioplex2_hygeo_bioid_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeo_bioid.featmat")

cilium_hygeo_df = pd.read_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_prey_pairs_pvalcorr_logchoose.feat")
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

cilium_hygeo_gb_df = cilium_hygeo_df.groupby(cilium_hygeo_df.index).first()
cilium_hygeo_avgspec2_gb_df = cilium_hygeo_avgspec2_df.groupby(cilium_hygeo_avgspec2_df.index).first()
cilium_hygeo_avgspec4_gb_df = cilium_hygeo_avgspec4_df.groupby(cilium_hygeo_avgspec4_df.index).first()

orig9k_bioplex2_hygeo_bioid_hygeo_df = orig9k_bioplex2_hygeo_bioid_df.join(cilium_hygeo_gb_df, how="outer", rsuffix="_cilium_hygeo")
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_df = orig9k_bioplex2_hygeo_bioid_hygeo_df.join(cilium_hygeo_avgspec2_gb_df, how="outer", rsuffix="_cilium_hygeo_avgspec2")
orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df = orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_df.join(cilium_hygeo_avgspec4_gb_df, how="outer", rsuffix="_cilium_hygeo_avgspec4")

orig9k_bioplex2_hygeo_bioid_hygeo_avgspec2_4_df.to_csv("/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/orig9k_bioplex2_hygeo_bioid_hygeo.featmat")



