# coding: utf-8
import pandas as pd
import numpy as np
import itertools as it
import scipy
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats

cilium_bioid_fname = "/stor/work/Marcotte/project/kdrew/data/cilium_bioid/P27_Prohits_viz_SAINT_input.txt"
cilium_bioid_df = pd.read_table(cilium_bioid_fname)

cilium_bioid_intatt_fname = "/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_suppTable1_Interactor_attributes_1-s2.0-S009286741501421X-mmc2-1.csv"
cilium_bioid_intatt_df = pd.read_csv(cilium_bioid_intatt_fname)
cilium_bioid_intatt_df = cilium_bioid_intatt_df.set_index("name")  

cilium_bioid_join_intatt_df = cilium_bioid_df.join(cilium_bioid_intatt_df, on="PreyGene", how="left")

#human_uniprot_df = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_20171218.tab")
#human_uniprot_df = human_uniprot_df.set_index(human_uniprot_df.columns[-1])
#cilium_bioid_uniprot_df = cilium_bioid_df.join(human_uniprot_df, on="PreyGene", how="left")

human_uniprot_df = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_reviewed_20171218.tab")
human_uniprot_df = human_uniprot_df.set_index(human_uniprot_df.columns[-1])
cilium_bioid_uniprot_df = cilium_bioid_df.join(human_uniprot_df, on="PreyGene", how="left")

preygenes = cilium_bioid_uniprot_df['PreyGene'].values
uniprot_genename_dict = human_uniprot_df['Gene names'].to_dict()

preygenes_genename = dict()
for g in preygenes:
    preygenes_genename[g] = None
    for k in uniprot_genename_dict:
        try:        
            if g in uniprot_genename_dict[k].split():
                    preygenes_genename[g] = k
                    break
        except AttributeError:
            continue

preygenes_genename_df = pd.DataFrame(preygenes_genename, index=['matched_genename']).T

preygenes_genename_nonan_df = preygenes_genename_df.query("matched_genename == matched_genename")

cilium_bioid_matched_df = cilium_bioid_df.join(preygenes_genename_nonan_df, on="PreyGene", how="left")
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="left")
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="left")

cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="inner")

cilium_bioid_matched_nonan_df = cilium_bioid_matched_df.query("matched_genename == matched_genename")
cilium_bioid_uniprot_df = cilium_bioid_matched_nonan_df.join(human_uniprot_df, on="matched_genename", how="left")

cilium_bioid_uniprot_nonan_df = cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] == cilium_bioid_uniprot_df['Cross-reference (GeneID)']]

cilium_bioid_uniprot_nonan_df['geneid'] = [x.split(';')[0] for x in cilium_bioid_uniprot_nonan_df['Cross-reference (GeneID)'].values]
cilium_bioid_uniprot_nonan_df.to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_uniprot_geneid.featmat", index=False)


cilium_bioid_uniprot_nonan_df['bait_genename'] = [x.split('_')[0] for x in cilium_bioid_uniprot_nonan_df['Bait'].values]

cilium_bioid_uniprot_baitgenename_df = cilium_bioid_uniprot_nonan_df.join(preygenes_genename_nonan_df, on="bait_genename", rsuffix='_bait', how="left")

bait_genename = dict()
for g in cilium_bioid_uniprot_nonan_df['bait_genename'].unique():
    bait_genename[g] = None
    for k in uniprot_genename_dict:
        try:        
            if g in uniprot_genename_dict[k].split():
                    bait_genename[g] = k
                    break
        except AttributeError:
            continue

bait_genename_df = pd.DataFrame(bait_genename, index=["matched_bait_genename"]).T

cilium_bioid_uniprot_bait_genename_df = cilium_bioid_uniprot_nonan_df.join(bait_genename_df, on="bait_genename", how="left")
cilium_bioid_uniprot_bait_df = cilium_bioid_uniprot_bait_genename_df.join(human_uniprot_df, on="matched_bait_genename", how="left", rsuffix="_bait")

cilium_bioid_uniprot_bait_df['bait_geneid'] = [x.split(';')[0] for x in cilium_bioid_uniprot_bait_df['Cross-reference (GeneID)_bait'].values]

cilium_bioid_uniprot_bait_df.to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_uniprot_geneid.featmat", index=False)


cilium_bioid_uniprot_bait_df['condition'] = ['_'.join(x.split('_')[1:]) for x in cilium_bioid_uniprot_bait_df.Bait.values]
cilium_bioid_featmat_df = cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']]

cilium_bioid_featmat_df['geneid_pairs_strsort'] = [str(sorted( [ str(x[0]), str(x[1]) ] )) for x in cilium_bioid_featmat_df[['bait_geneid','geneid']].values]

cilium_bioid_featmat_df = cilium_bioid_featmat_df.set_index("geneid_pairs_strsort")

cilium_bioid_featmat_df.query("condition == 'ciliated'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_ciliated.featmat")
cilium_bioid_featmat_df.query("condition == 'non_ciliated'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_non_ciliated.featmat")
cilium_bioid_featmat_df.query("condition == 'FLAG_IP'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_FLAG_IP.featmat")
cilium_bioid_featmat_df.to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_All.featmat")
cilium_bioid_featmat_df.query("condition == ''").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_none.featmat")

cilium_bioid_hypergeo_input_df = cilium_bioid_featmat_df.reset_index()[['Bait','geneid','AvgSpec','condition']]
cilium_bioid_hypergeo_input_df['dataset'] = [x+"_bioid" for x in cilium_bioid_hypergeo_input_df['condition'].values]

from scipy import stats
cilium_bioid_hypergeo_bioidonly_df = cilium_bioid_hypergeo_input_df.query("condition == 'ciliated' or condition == 'non_ciliated'")
pzscores = (stats.rankdata(cilium_bioid_hypergeo_bioidonly_df['AvgSpec'])-1) / (len(cilium_bioid_hypergeo_bioidonly_df['AvgSpec'])-1) * 100
cilium_bioid_hypergeo_bioidonly_df['abundance'] = pzscores
cilium_bioid_hypergeo_bioidonly_df.columns = ['experiment_id','geneid','AvgSpec','condition','dataset','abundance']
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input.featmat")

cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 2")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input_AvgSpec2.featmat")
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 4")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input_AvgSpec4.featmat")

