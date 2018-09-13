# coding: utf-8
import pandas as pd
orig9k_fname = "/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_feature_revisitMerge_pairsOnly_preyMerge2_heinCollapseMerge_pairsOnly_preyMerge2.txt"
orig9k_df = pd.read_csv(orig9k_fname)
full_featmat_fname = "/stor/work/Marcotte/project/kdrew/data/protein_complex_maps/complex_map2/orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean.featmat"
full_featmat_df = pd.read_csv(full_featmat_fname)
orig9k_df.head()
full_featmat_df.head()
orig9k_df.tail()
full_featmat_df.tail()
full_featmat_df.head(2)
full_featmat_df.head(2).values
orig9k_df.head(2)
orig9k_df.query("gene_id1_y == gene_id1_y").head(2)
orig9k_df.query("gene_id1_y == gene_id1_y")['enrichment']
orig9k_df.query("enrichment == enrichment").head()
orig9k_df.query("enrichment == enrichment").shape
orig9k_df.query("enrichment == enrichment and gene_id1_y == gene_id1_y").shape
orig9k_df.query("enrichment == enrichment and gene_id1_y == gene_id1_y and gene_id2_y == gene_id2_y").shape
orig9k_enrichment_df = orig9k_df.query("enrichment == enrichment and gene_id1_y == gene_id1_y and gene_id2_y == gene_id2_y")
orig9k_enrichment_df = orig9k_enrichment_df.set_index("frozenset_geneids")
full_featmat_df.head(2)
full_featmat_df.columns
list(full_featmat_df.columns)
full_featmat_df.head(2)
full_featmat_df.head(2).values
full_featmat_df.head(2)
pairs = full_featmat_df[['id1','id2']].values
pairs[:10]
pairs_strsort = [str(sorted([str(int(x[0])),str(int(x[1]))])) for x in pairs]
pairs_strsort = [str(sorted([str(x[0]),str(x[1])])) for x in pairs]
pairs_strsort[:10]
pairs_strsort[-10:]
pairs_strsort[1000:1010]
full_featmat_df['pairs_strsort'] = pairs_strsort
full_featmat_df = full_featmat_df.set_index('pairs_strsort')
full_featmat_df.head(2)
orig9k_enrichment_df.head(2)
enrichment_pairs = orig9k_enrichment_df[['gene_id1_str_y','gene_id2_str_y']].values
orig9k_enrichment_df['pairs_strsort'] = [str(sorted([str(int(x[0])),str(int(x[1]))])) for x in enrichment_pairs]
orig9k_enrichment_df.head(2)
orig9k_enrichment_df = orig9k_enrichment_df.set_index('pairs_strsort')
orig9k_enrichment_df.head(2)
orig9k_enrichment_df['enrichment'].head()
full_featmat_wEnrichment_df = full_featmat_df.join(orig9k_enrichment_df['enrichment'], how='left')
full_featmat_wEnrichment_df.head(2)
full_featmat_wEnrichment_df.query("enrichment == enrichment").shape
full_featmat_wEnrichment_df.to_csv("/project/kdrew/data/protein_complex_maps/complex_map2.5/orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_trimCols_groupbyMean_wHeinEnrichment.featmat",index=False) 
get_ipython().magic(u'save add_hein_enrichment_featmat.py 1:51')
