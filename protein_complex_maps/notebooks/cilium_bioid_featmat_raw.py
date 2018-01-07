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
cilium_bioid_df.head()
cilium_bioid_intatt_df = "/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_suppTable1_Interactor_attributes_1-s2.0-S009286741501421X-mmc2-1.csv"
cilium_bioid_intatt_fname = "/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_suppTable1_Interactor_attributes_1-s2.0-S009286741501421X-mmc2-1.csv"
cilium_bioid_intatt_df = pd.read_csv(cilium_bioid_intatt_fname)
cilium_bioid_intatt_df.head()
cilium_bioid_intatt_df.head(1)
cilium_bioid_df.head()
cilium_bioid_intatt_df.head(2)
cilium_bioid_int_df = cilium_bioid_intatt_df.set_index("name")
cilium_bioid_intatt_df = cilium_bioid_intatt_df.set_index("name")  
cilium_bioid_join_intatt_df = cilium_bioid_df.join(cilium_bioid_intatt_df, on="PreyGene", how="left")
cilium_bioid_join_intatt_df.head(1)
cilium_bioid_join_intatt_df[['ENTREZ gene ID']].head()
cilium_bioid_join_intatt_df[['PreyGene','ENTREZ gene ID']].head()
cilium_bioid_join_intatt_df[['ENTREZ gene ID']].head()
cilium_bioid_join_intatt_df[['ENTREZ gene ID']].isnan()
cilium_bioid_join_intatt_df[['ENTREZ gene ID']].values
cilium_bioid_join_intatt_df[['PreyGene','ENTREZ gene ID']].tail()
human_uniprot_df = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_20171218.tab")
human_uniprot_df.head()
human_uniprot_df = human_uniprot_df.set_index("Gene names (primary )")
human_uniprot_df.columns
human_uniprot_df.columns[-1]
human_uniprot_df = human_uniprot_df.set_index(human_uniprot_df.columns[-1])
cilium_bioid_uniprot_df = cilium_bioid_df.join(human_uniprot_df, on="PreyGene", how="left")
cilium_bioid_uniprot_df.head()
human_uniprot_df = pd.read_table("/project/kdrew/data/uniprot/uniprot-proteome_hsapiens_reviewed_20171218.tab")
human_uniprot_df = human_uniprot_df.set_index(human_uniprot_df.columns[-1])
cilium_bioid_uniprot_df = cilium_bioid_df.join(human_uniprot_df, on="PreyGene", how="left")
cilium_bioid_uniprot_df.head()
cilium_bioid_uniprot_df.query("Cross-reference (GeneID) == Cross-reference (GeneID)")
cilium_bioid_uniprot_df.query("'Cross-reference (GeneID)' != 'Cross-reference (GeneID)'")
cilium_bioid_uniprot_df.head()
cilium_bioid_uniprot_df.tail()
cilium_bioid_uniprot_df.tail(1)
cilium_bioid_uniprot_df['Bait']
cilium_bioid_uniprot_df.tail(1)
cilium_bioid_uniprot_df.tail(2)
cilium_bioid_uniprot_df['Cross-reference (GeneID)']
cilium_bioid_uniprot_df['Cross-reference (GeneID)'] == cilium_bioid_uniprot_df['Cross-reference (GeneID)']
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
get_ipython().magic(u'save cilium_bioid_featmat.py 1:51')
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']].head()
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]['PreyGene']
len(cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]['PreyGene'])
nan_preygenes = cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]['PreyGene']
human_uniprot_df.head()
human_uniprot_df['Gene names']
human_uniprot_df['Gene names'].values
genename_lists = human_uniprot_df['Gene names'].values
nan_preygenes
nan_preygenes.values
nan_preygenes.values[0]
nan_preygenes.values[0] in "xxx C19orf43 yyy"
[i if 'C19orf43' in genenames for i, genenames in enumerate(genename_lists)]
[i if 'C19orf43' in genenames else -1 for i, genenames in enumerate(genename_lists)]
[i for i, genenames in enumerate(genename_lists)]
[nan_preygenes[0] in genenames for genenames in genename_lists]
nan_preygenes
nan_preygenes.values[0]
[nan_preygenes.values[0] in genenames for genenames in genename_lists]
nan_preygenes
cilium_bioid_uniprot_df['PreyGene']
cilium_bioid_uniprot_df['PreyGene'].apply(lambda x: human_uniprot_df['
cilium_bioid_uniprot_df.query("PreyGene in ' YARS '")
cilium_bioid_uniprot_df.query("PreyGene in 'YARS'")
cilium_bioid_uniprot_df.query("PreyGene in ' YARS'")
cilium_bioid_uniprot_df.query("PreyGene == 'YARS'")
cilium_bioid_uniprot_df.query("' YARS '.str.contains(PreyGene)")
cilium_bioid_uniprot_df['PreyGene']
preygenes = cilium_bioid_uniprot_df['PreyGene'].values() 
preygenes = cilium_bioid_uniprot_df['PreyGene'].values
human_uniprot_df.head()
human_uniprot_df[['Gene names  (primary )']]
human_uniprot_df.head()
human_uniprot_df[['index','Gene names']]
human_uniprot_df['Gene names']
list(human_uniprot_df['Gene names'])
list(human_uniprot_df['Gene names'])[0]
human_uniprot_df['Gene names'].todict()
human_uniprot_df['Gene names'].to_dict()
uniprot_genename_dict = human_uniprot_df['Gene names'].to_dict()
preygenes
[i[0] for i in uniprot_genename_dict.items()]
[i[0] for i in uniprot_genename_dict.items() if g in i[1] for g in preygenes]
[i[0] for i in uniprot_genename_dict.items() for g in preygenes if g in i[1]]
preygenes
[x for x in preygenes]
[i[0] for i in uniprot_genename_dict.items() for g in preygenes if g in i[1]]
[i[0] if g in i[1] for g in preygenes for i in uniprot_genename_dict.items()]
[i[0] for g in preygenes if g in i[1] for i in uniprot_genename_dict.items()]
preygenes_genename = []
for g in preygenes:
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k]:
			preygenes_genename.append(k)
			break
for k in uniprot_genename_dict:
	print k
for g in genenames:
	print g
genenames
for g in preygenes:
	print g
for g in preygenes:
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k]:
			preygenes_genename.append(k)
			break
uniprot_genename_dict['YARS']
'YARS' in uniprot_genename_dict['YARS'] 
for g in preygenes:
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k]:
			preygenes_genename.append(k)
			break
preygenes_genename = dict()
for g in preygenes:
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k]:
			preygenes_genename[g] = k
			break
for g in preygenes:
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k]:
			preygenes_genename[g] = k
			break
get_ipython().magic(u'edit')
for g in preygenes:
	for k in uniprot_genename_dict: preygenes_genename
		if g in uniprot_genename_dict[k]:
			preygenes_genename.append(k)
			break
preygenes_genename
for g in preygenes:
	preygenes_genename[g] = None
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k].split()
			preygenes_genename[g] = k
			break
for g in preygenes:
	preygenes_genename[g] = None
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k].split():
			preygenes_genename[g] = k
			break
for g in preygenes:
	preygenes_genename[g] = None
	for k in uniprot_genename_dict:
		if g in uniprot_genename_dict[k].split():
			preygenes_genename[g] = k
			break
get_ipython().magic(u'save cilium_bioid_featmat.py 1:120')
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
preygenes_genename.keys()[0]
preygenes_genename[preygenes_genename.keys()[0]]
preygenes_genename[preygenes_genename.keys()[1]]
print preygenes_genename.keys()[1]; preygenes_genename[preygenes_genename.keys()[1]]
i = 2; print preygenes_genename.keys()[i]; preygenes_genename[preygenes_genename.keys()[i]]
i = 3; print preygenes_genename.keys()[i]; preygenes_genename[preygenes_genename.keys()[i]]
preygene_genename['RPS10-NUDT3']
preygenes_genename['RPS10-NUDT3']
preygenes_genename['RPS10']
preygenes_genename['NUDT3']
nan_preygenes[0]
nan_preygenes.head()
preygenes_genename['C19orf43']
preygenes_genename['C22orf28']
preygenes_genename['BOLA2']
preygenes_genename['ATP5J2-PTCD1']
preygenes_genename['56160532']
nan_preygenes.head(10)
preygenes_genename['LOC645870']
preygenes_genename['ATP6']
preygenes_genename['NHP2L1']
preygenes_genename['MESDC2']
preygenes_genename['FAM158A']
preygenes_genename_df = pd.DataFrame(preygenes_genename).T
preygenes_genename_df = pd.DataFrame(preygenes_genename)
preygenes_genename_df = pd.DataFrame.from_dict(preygenes_genename)
preygenes_genename_df = pd.DataFrame().from_dict(preygenes_genename)
preygenes_genename_df = pd.DataFrame().from_dict(preygenes_genename, index=[0])
preygenes_genename_df = pd.DataFrame(preygenes_genename, index=[0])
preygenes_genename_df.head()
preygenes_genename_df = pd.DataFrame(preygenes_genename, index=[0]).T
preygenes_genename_df.head()
preygenes_genename_df.head(10)
preygenes_genename_df = pd.DataFrame(preygenes_genename, index=['matched_genename']).T
preygenes_genename_df.head(10)
cilium_bioid_uniprot_df.head()
cilium_bioid_df.head()
cilium_bioid_matched_df = cilium_bioid_df.join(preygenes_genename_df, on="PreyGene", how="left")
cilium_bioid_matched_df.head()
cilium_bioid_matched_df.query("PreyGene == 'FAM158A')
cilium_bioid_matched_df.query("PreyGene == 'FAM158A'")
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="left")
cilum_bioid_uniprot_df.head()
cilium_bioid_uniprot_df.head()
cilium_bioid_uniprot_df.query("PreyGene == 'FAM158A'")
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']].head()
preygenes_genename_df.head()
preygenes_genename_df.query("matched_genename != None")
preygenes_genename_df.query("matched_genename != matched_genename")
preygenes_genename_df.query("matched_genename == matched_genename")
preygenes_genename_nonan_df = preygenes_genename_df.query("matched_genename == matched_genename")
cilium_bioid_matched_df = cilium_bioid_df.join(preygenes_genename_nonan_df, on="PreyGene", how="left")
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="left")
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']].head()
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="left")
human_uniprot_df.head()
human_uniprot_df.index == human_uniprot_df.index]
human_uniprot_df[human_uniprot_df.index != human_uniprot_df.index]
human_uniprot_df[human_uniprot_df.index == human_uniprot_df.index]
get_ipython().magic(u'reset ')
cilium_bioid_uniprot_df.head()
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']].head()
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']].head(10)
cilium_bioid_df.head()
cilium_bioid_df.query("Prey == 'ATP5J2-PTCD1'")
cilium_bioid_df.query("PreyGene == 'ATP5J2-PTCD1'")
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
cilium_bioid_uniprot_df
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
cilium_bioid_df
cilium_bioid_matched_df
cilium_bioid_uniprot_df = cilium_bioid_matched_df.join(human_uniprot_df, on="matched_genename", how="inner")
cilium_bioid_uniprot_df
cilium_bioid_matched_df
cilium_bioid_matched_df.query("matched_genename == matched_genename")
cilium_bioid_matched_df.query("matched_genename != matched_genename")
cilium_bioid_matched_nonan_df = cilium_bioid_matched_df.query("matched_genename == matched_genename")
cilium_bioid_matched_nonan_df
cilium_bioid_uniprot_df = cilium_bioid_matched_nonan_df.join(human_uniprot_df, on="matched_genename", how="left")
cilium_bioid_uniprot_df
cilium_bioid_uniprot_df['Cross-reference (GeneID)']
cilium_bioid_uniprot_df['Cross-reference (GeneID)'].values
[x.split(';')[0] for x in cilium_bioid_uniprot_df['Cross-reference (GeneID)'].values]
cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] != cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
cilium_bioid_uniprot_nonan_df = cilium_bioid_uniprot_df[cilium_bioid_uniprot_df['Cross-reference (GeneID)'] == cilium_bioid_uniprot_df['Cross-reference (GeneID)']]
[x.split(';')[0] for x in cilium_bioid_uniprot_nonan_df['Cross-reference (GeneID)'].values]
cilium_bioid_uniprot_nonan_df['geneid'] = [x.split(';')[0] for x in cilium_bioid_uniprot_nonan_df['Cross-reference (GeneID)'].values]
cilium_bioid_uniprot_nonan_df['geneid'] = [x.split(';')[0] for x in cilium_bioid_uniprot_nonan_df['Cross-reference (GeneID)'].values]
cilium_bioid_uniprot_nonan_df.head()
cilium_bioid_uniprot_nonan_df.to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_uniprot_geneid.featmat", index=False)
cilium_bioid_uniprot_nonan_df.head()
cilium_bioid_uniprot_nonan_df['Bait']
cilium_bioid_uniprot_nonan_df['Bait'].values
[x.split('_')[0] for x in cilium_bioid_uniprot_nonan_df['Bait'].values]
cilium_bioid_uniprot_non_df['bait_genename'] = [x.split('_')[0] for x in cilium_bioid_uniprot_nonan_df['Bait'].values]
cilium_bioid_uniprot_nonan_df['bait_genename'] = [x.split('_')[0] for x in cilium_bioid_uniprot_nonan_df['Bait'].values]
cilium_bioid_uniprot_nonan_df.head()
preygenes_genename_df.head()
preygenes_genename_nonan_df.head()
cilium_bioid_uniprot_baitgenename_df = cilium_bioid_uniprot_nonan_df.join(preygenes_genename, on="bait_genename", rsuffix='_bait', how="left")
cilium_bioid_uniprot_baitgenename_df = cilium_bioid_uniprot_nonan_df.join(preygenes_genename, on=["bait_genename"], rsuffix='_bait', how="left")
preygenes_genename.head()
cilium_bioid_uniprot_baitgenename_df = cilium_bioid_uniprot_nonan_df.join(preygenes_genename_nonan_df, on="bait_genename", rsuffix='_bait', how="left")
cilium_bioid_uniprot_baitgenename_df.head()
cilium_bioid_uniprot_baitgenename_df.query("matched_genename_bait != matched_genename_bait")
cilium_bioid_uniprot_nonan_df.head()
cilium_bioid_uniprot_nonan_df['bait_genename'].values
cilium_bioid_uniprot_nonan_df['bait_genename'].unique()
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
bait_genename
get_ipython().magic(u'save cilium_bioid_featmat_raw 1:236')
bait_genename
len(bait_genename)
bait_genename_df = pd.DataFrame(bait_genename, index=["matched_bait_genename"])
bait_genename_df
bait_genename_df = pd.DataFrame(bait_genename, index=["matched_bait_genename"]).T
bait_genename_df
cilium_bioid_uniprot_bait_genename_df = cilium_bioid_uniprot_nonan_df.join(bait_genename_df, on="bait_genename", how="left")
cilium_bioid_uniprot_bait_genename_df.head()
cilium_bioid_uniprot_bait_df = cilium_bioid_uniprot_bait_genename_df.join(human_uniprot_df, on="matched_bait_genename", how="left")
cilium_bioid_uniprot_bait_df = cilium_bioid_uniprot_bait_genename_df.join(human_uniprot_df, on="matched_bait_genename", how="left", rsuffix="_bait")
cilium_bioid_uniprot_bait_df.head()
cilium_bioid_uniprot_bait_df[cilium_bioid_uniprot_bait_df['Cross-reference (GeneID)_bait'] != cilium_bioid_uniprot_bait_df['Cross-reference (GeneID)_bait']]
cilium_bioid_uniprot_bait_df['bait_geneid'] = [x.split(';')[0] for x in cilium_bioid_uniprot_bait_df['Cross-reference (GeneID)_bait'].values]
cilium_bioid_uniprot_bait_df.head()
cilium_bioid_uniprot_bait_df.to_csv("cilium_bioid_uniprot_geneid.featmat", index=False)
cilium_bioid_uniprot_bait_df.to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_uniprot_geneid.featmat", index=False)
cilium_bioid_uniprot_bait_df.columns
cilium_bioid_uniprot_bait_df.Bait.unique()
[x.split('_')[1:] for x in cilium_bioid_uniprot_bait_df.Bait.values]
cilium_bioid_uniprot_bait_df['condition'] = [x.split('_')[1:] for x in cilium_bioid_uniprot_bait_df.Bait.values]
cilium_bioid_uniprot_bait_df.head()
cilium_bioid_uniprot_bait_df['condition'] = ['_'.join(x.split('_')[1:]) for x in cilium_bioid_uniprot_bait_df.Bait.values]
cilium_bioid_uniprot_bait_df.head()
cilium_bioid_uniprot_bait_df.columns
cilium_bioid_uniprot_bait_df.condition.unique()
cilium_bioid_uniprot_bait_df.query("AvgP != TopoAvgP")
cilium_bioid_uniprot_bait_df.query("AvgP != TopoAvgP")[['Bait','PreyGene','AvgP','TopoAvgP']]
cilium_bioid_uniprot_bait_df.query("MaxP != TopoMaxP")[['Bait','PreyGene','AvgP','TopoAvgP','MaxP','TopoMaxP']]
cilium_bioid_uniprot_bait_df.query("AvgP != TopoAvgP")[['Bait','PreyGene','AvgP','TopoAvgP','MaxP','TopoMaxP','SaintScore']]
cilium_bioid_uniprot_bait_df.query("AvgP != SaintScore")[['Bait','PreyGene','AvgP','TopoAvgP','MaxP','TopoMaxP','SaintScore']]
126+120
126.0/120
cilium_bioid_uniprot_bait_df.query("AvgP != SaintScore")[['Bait','PreyGene','AvgP','TopoAvgP','MaxP','TopoMaxP','SaintScore']]
cilium_bioid_uniprot_bait_df.columns
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].head()
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].head()
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].sort_values("AvgSpec").head()
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].sort_values("AvgSpec")
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].query("AvgSpec == 1")
cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']].head()
cilium_bioid_featmat_df = cilium_bioid_uniprot_bait_df[['Bait','PreyGene','AvgSpec','AvgP','MaxP','Fold_Change','BFDR','bait_geneid','geneid','condition']]
cilium_bioid_featmat_df.head()
[sorted(list(str(x[0]),str(x[1]))) for x in cilium_bioid_featmat_df[['bait_geneid','geneid']].values]
[sorted(list((str(x[0]),str(x[1])))) for x in cilium_bioid_featmat_df[['bait_geneid','geneid']].values]
[str(sorted( [ str(x[0]), str(x[1]) ] )) for x in cilium_bioid_featmat_df[['bait_geneid','geneid']].values]
cilium_bioid_featmat_df['geneid_pairs_strsort'] = [str(sorted( [ str(x[0]), str(x[1]) ] )) for x in cilium_bioid_featmat_df[['bait_geneid','geneid']].values]
cilium_bioid_featmat_df.head()
cilium_bioid_featmat_df = cilium_bioid_featmat_df.set_index("geneid_pairs_strsort")
cilium_bioid_featmat_df.head()
cilium_bioid_featmat_df.condition.unique()
cilium_bioid_featmat_df.query("condition == 'ciliated'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_ciliated.featmat")
cilium_bioid_featmat_df.query("condition == 'non_ciliated'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_non_ciliated.featmat")
cilium_bioid_featmat_df.query("condition == 'FLAG_IP'").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_FLAG_IP.featmat")
cilium_bioid_featmat_df.to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_All.featmat")
cilium_bioid_featmat_df.query("condition == ''").to_csv("/stor/work/Marcotte/project/kdrew/data/cilium_bioid/cilium_bioid_none.featmat")
get_ipython().magic(u'save cilium_bioid_featmat_raw 1:291')
cilium_bioid_uniprot_geneid.head()
cilium_bioid_featmat_df.head()
cilium_bioid_hypergeo_intput_df = cilium_bioid_featmat_df[['Bait','geneid','AvgSpec']]
cilium_bioid_hypergeo_input_df = cilium_bioid_featmat_df[['Bait','geneid','AvgSpec']]
cilium_biolid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df.reset_index()
cilium_bioid_hypergeo_input_df = cilium_bioid_featmat_df.reset_index()[['Bait','geneid','AvgSpec']]
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df.columns = ['experiment_id','geneid','abundance']
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df['dataset'] = "cilium_bioid"
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_featmat_df.head()
cilium_bioid_hypergeo_input_df = cilium_bioid_featmat_df[['Bait','geneid','AvgSpec','condition']]
cilium_bioid_hypergeo_input_df.columns = ['experiment_id','geneid','abundance','condition']
cilium_bioid_hypergeo_input_df[cilium_bioid_hypergeo_input_df['condition'] == 'ciliated'].head()
cilium_bioid_hypergeo_input_df[cilium_bioid_hypergeo_input_df['condition'] == 'ciliated']['dataset'] = "ciliated_bioid"
cilium_bioid_hypergeo_input_df[cilium_bioid_hypergeo_input_df['condition'] == 'non_ciliated']['dataset'] = "non_ciliated_bioid"
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df = cilium_bioid_featmat_df.reset_index()[['Bait','geneid','AvgSpec','condition']]
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df.query("condition == ciliated")['data_set'] = 'ciliated_bioid'
cilium_bioid_hypergeo_input_df.query("condition == 'ciliated'")['data_set'] = 'ciliated_bioid'
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df['dataset'] = ''
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df.query("condition == 'ciliated'")['dataset'] = 'ciliated_bioid'
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_input_df['dataset'] = [x+"_bioid" for x in cilium_bioid_hypergeo_input_df['condition'].values]
cilium_bioid_hypergeo_input_df.head()
from scipy import stats
cilium_bioid_hypergeo_input_df.head()
cilium_bioid_hypergeo_bioidonly_df = cilium_bioid_hypergeo_input_df.query("condition == 'ciliated' or condition == 'non_ciliated')
cilium_bioid_hypergeo_bioidonly_df = cilium_bioid_hypergeo_input_df.query("condition == 'ciliated' or condition == 'non_ciliated'")
cilium_bioid_hypergeo_bioidonly_df.head()
len(cilium_bioid_hypergeo_bioidonly_df)
len(cilium_bioid_hypergeo_input_df)
pzscores = (stats.rankdata(cilium_bioid_hypergeo_bioidonly_df['AvgSpec'])-1) / (len(cilium_bioid_hypergeo_bioidonly_df['AvgSpec'])-1) * 100
pzscores[:10]
cilium_bioid_hypergeo_input_df['abundance'] = pzscores
cilium_bioid_hypergeo_bioidonly_df['abundance'] = pzscores
cilium_bioid_hypergeo_bioidonly_df.head()
cilium_bioid_hypergeo_bioidonly_df.columns = ['experiment_id','geneid','abundance','dataset']
cilium_bioid_hypergeo_bioidonly_df.columns = ['experiment_id','geneid','AvgSpec','condition','dataset','abundance']
cilium_bioid_hypergeo_bioidonly_df.head()
cilium_bioid_hypergeo_bioidonly_df['experiment_id','geneid','abundance','dataset'].head()
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']]
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']].head()
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']].dataset.unique()
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']].head()
cilium_bioid_hypergeo_bioidonly_df[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input.featmat")
cilium_bioid_hypergeo_bioidonly_df.sorted_values('abundance').head()
cilium_bioid_hypergeo_bioidonly_df.sort_values('abundance').head()
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 2")[['experiment_id','geneid','abundance','dataset']].to_csv("/poject/kdrew/data/cilium_bioid/cilium_bioid_fhypergeo_input.featmat")
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 2")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_fhypergeo_inpuhypergeo_input_AvgSpec2.featmat")
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 2")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input_AvgSpec2.featmat")
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 4")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/"cilium_bioid_hypergeo_input_AvgSpec4.featmat")
cilium_bioid_hypergeo_bioidonly_df.query("AvgSpec >= 4")[['experiment_id','geneid','abundance','dataset']].to_csv("/project/kdrew/data/cilium_bioid/cilium_bioid_hypergeo_input_AvgSpec4.featmat")
