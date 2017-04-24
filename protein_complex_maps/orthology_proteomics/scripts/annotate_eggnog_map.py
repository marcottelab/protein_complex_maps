from __future__ import print_function
import pandas as pd
import sys
import argparse


def make_labels(eggnog, outfile):

   #output from get_eggion_ids.py
   egg=pd.read_csv(eggnog, index_col=False, sep="\t")
   
   egg= egg[['GroupID', 'ProteinID', 'Annotation']]

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
   all_prots = egg.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
   all_prots.columns= ['GroupID', 'Eggnog_annotation', 'AllMembers']


   #print all_prots
   one_prot = egg.groupby(['GroupID',  'Annotation']).head(1)
   one_prot = one_prot[['GroupID', 'ProteinID']]
   

 
   all_prots = all_prots.set_index(["GroupID"])

   
   one_prot = one_prot.set_index(["GroupID"])

   #print all_prots
   #print one_prot


   full_table = all_prots.join(one_prot, how = "left")
   split_table = full_table['ProteinID'].str.split("|", return_type='frame')
   print(split_table)
   #get accession number
   split_table=split_table[[1]]

   #if key == 'acc' or key == 'entry':
   #      split_table = pd.DataFrame(df.row.str.split('|',1).tolist(),
   #                                  columns = ['cat','acc', 'entry'])



   final = full_table.join(split_table, how="left")
   final = final[['Eggnog_annotation', 'AllMembers', 1]]
   final.columns= ['Eggnog_annotation', 'AllMembers', 'Accession'] 
   
   #print full_table

   final.to_csv(outfile, sep="\t")



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

    
parser = argparse.ArgumentParser(description='Short sample app')

#parser.add_argument('identified_eggion', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
parser.add_argument('outfile', action="store", type=str)
inputs = parser.parse_args()

make_labels(inputs.orthology_file, inputs.outfile)






    
