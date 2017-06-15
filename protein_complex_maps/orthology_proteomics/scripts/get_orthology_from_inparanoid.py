from __future__ import print_function
import pandas as pd
import sys
import argparse
'''
1       4369    C.reinhardtii   1.000   A8JHU7  100%
1       4369    H.sapiens       1.000   Q8IVF4  100%
2       4074    C.reinhardtii   1.000   A8IUF0  100%
2       4074    H.sapiens       1.000   Q9P225  100%
3       4008    C.reinhardtii   1.000   A8HPL0  100%
3       4008    H.sapiens       1.000   Q6P2Q9  100%
4       3901    C.reinhardtii   1.000   A8J063  100%
4       3901    H.sapiens       1.000   Q9P2D7  100%
5       3660    C.reinhardtii   1.000   A8I5C3  95%
5       3660    H.sapiens       1.000   Q8WXX0  100%
'''

def get_pairwise_from_sql(inparanoid_orthology):
   
    inparanoid = pd.read_csv(inparanoid_orthology, index_col=False, header=None, sep='\t')
    
    inparanoid.columns = ['id1', 'id2', 'spec', 'score', 'acc', 'percent']
    ip_grouped = inparanoid.groupby(['id1','spec']).head(1)
    
    print(ip_grouped)
    
    selection = ip_grouped[['id1', 'id2', 'spec', 'acc']]
    selection = selection.reset_index() 
    odd_df = selection[selection.index % 2 != 0]
    print(odd_df[odd_df.id1==16])
    odd_df = odd_df.set_index(['id1'])
    #odd_df['pair1'] = odd_df[['spec', 'acc']].apply(lambda x: "|".join(x), axis=1)
    #odd_df = odd_df.drop(['spec', 'acc'], 1)
  


    #print(odd_df)
    even_df = selection[selection.index % 2  == 0]
    print(even_df[even_df.id1==16])
    even_df = even_df.set_index(['id1'])
    df = odd_df.join(even_df, how="outer", rsuffix='_right')
    final_df = df[['spec', 'acc', 'spec_right', 'acc_right']]
    final_df.columns = ['spec_1', 'Entry', 'spec_2', 'acc_2']

    print(final_df)
    outfile = "conversiontbl_" + inparanoid_orthology
    final_df.to_csv(outfile, index=False)

    #dfall['ID1'] = [list(x)[0] for x in dfall['pairset'].values]
 
parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('inparanoid_orthology', action="store", type=str)
inputs = parser.parse_args()


get_pairwise_from_sql(inputs.inparanoid_orthology)






'''
def make_labels(eggnog, profile, outfile):

   #output from get_elution_ids.py
   egg=pd.read_csv(eggnog, index_col=False, sep="\t")
   
   egg= egg[['GroupID', 'AllMembers']]

   egg = egg.set_index(['GroupID'])

   prof = pd.read_csv(profile, index_col=False, sep=",")
   prof = prof.set_index(['GroupID'])

   annotated = egg.join(prof, how="right")
   print annotated 



   #changing from long to wide format for eggion profiles
#   wide = egg.pivot(index='GroupID', columns = 'FractionID', values='Total_SpecCounts')  
   #print wide

#   wide = wide.fillna(0)


#   raw_outfile = identified_eggion.replace("_eggion_", "_raw_wide_eggion_")
#   wide.to_csv(raw_outfile)

 #  wide['Total'] = wide.sum(axis=1)


   #saving and importing csv fixes header columns from pivot format
 #  wide2 = pd.read_csv(raw_outfile)
 #  wide2 = wide2.set_index(['GroupID'])

   #Pull annotations from orthology file (eggnog_output)
 #  annot = pd.read_csv(orthology_file, sep="\t")
 #  annot = annot[['GroupID', 'ProteinID', 'Annotation']]
   #annot = annot.set_index(['GroupID']) 

   #Consolidate proteins from a group into one row
   #all_prots = egg.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
   #all_prots.columns= ['GroupID', 'Eggnog_annotation', 'AllMembers']


   #print all_prots
   #one_prot = egg.groupby(['GroupID',  'Annotation']).head(1)
   #one_prot = one_prot[['GroupID', 'ProteinID']]
   

 
   #all_prots = all_prots.set_index(["GroupID"])

   
   #one_prot = one_prot.set_index(["GroupID"])

   #print all_prots
   #print one_prot


   #full_table = all_prots.join(one_prot, how = "left")
   #split_table = full_table['ProteinID'].str.split("|", return_type='frame')
   #split_table=split_table[[2]]

   #final = full_table.join(split_table, how="left")
   #final = final[['Eggnog_annotation', 'AllMembers', 2]]
   #final.columns= ['Eggnog_annotation', 'AllMembers', 'Entry']

 


   #print full_table
   #outfile = eggnog.replace("_orthology", "_map_annotations")


   #print outfile
   annotated.to_csv(outfile, sep=",")



   #Get annotations and 
   #ungrouped_alt_wide = annot.join(wide, how = "right")

   #ungrouped_alt_wide = ungrouped_alt_wide.reset_index()
   #http://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
   #alt_wide_labels = ungrouped_alt_wide.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
#   print grouped_egg

#   alt_wide_labels = alt_wide_labels.set_index(["GroupID"])
#   alt_wide = alt_wide_labels.join(wide, how = "inner")
#   print alt_wide
   #print alt_wide

 #  alt_outfile = identified_eggion.replace("_eggion_", "_alt_wide_eggion_")
 #  alt_wide.to_csv(alt_outfile)
'''
   
    
