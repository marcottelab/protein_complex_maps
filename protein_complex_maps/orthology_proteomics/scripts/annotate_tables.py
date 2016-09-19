import pandas as pd
import sys
import argparse


def make_labels(annotation, profile, outfile, join_on):

   #output from get_elution_ids.py
   annot=pd.read_csv(annotation, index_col=False, sep="\t")
  

   annot = annot.set_index([join_on])

   prof = pd.read_csv(profile, index_col=False, sep="\t")
   prof = prof.set_index([join_on])

   annotated = prof.join(annot, how="left")
   print annotated 



   #changing from long to wide format for annotion profiles
#   wide = annot.pivot(index='GroupID', columns = 'FractionID', values='Total_SpecCounts')  
   #print wide

#   wide = wide.fillna(0)


#   raw_outfile = identified_annotion.replace("_annotion_", "_raw_wide_annotion_")
#   wide.to_csv(raw_outfile)

 #  wide['Total'] = wide.sum(axis=1)


   #saving and importing csv fixes header columns from pivot format
 #  wide2 = pd.read_csv(raw_outfile)
 #  wide2 = wide2.set_index(['GroupID'])

   #Pull annotations from annotation file (annotnog_output)
 #  annot = pd.read_csv(annotation_file, sep="\t")
 #  annot = annot[['GroupID', 'ProteinID', 'Annotation']]
   #annot = annot.set_index(['GroupID']) 

   #Consolidate proteins from a group into one row
   #all_prots = annot.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
   #all_prots.columns= ['GroupID', 'Eggnog_annotation', 'AllMembers']


   #print all_prots
   #one_prot = annot.groupby(['GroupID',  'Annotation']).head(1)
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
   #outfile = annotnog.replace("_annotation", "_map_annotations")


   #print outfile
   annotated.to_csv(outfile, sep=",")



   #Get annotations and 
   #ungrouped_alt_wide = annot.join(wide, how = "right")

   #ungrouped_alt_wide = ungrouped_alt_wide.reset_index()
   #http://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
   #alt_wide_labels = ungrouped_alt_wide.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
#   print grouped_annot

#   alt_wide_labels = alt_wide_labels.set_index(["GroupID"])
#   alt_wide = alt_wide_labels.join(wide, how = "inner")
#   print alt_wide
   #print alt_wide

 #  alt_outfile = identified_annotion.replace("_annotion_", "_alt_wide_annotion_")
 #  alt_wide.to_csv(alt_outfile)

    
parser = argparse.ArgumentParser(description='Short sample app')

#parser.add_argument('identified_annotion', action="store", type=str)
parser.add_argument('annotation_file', action="store", type=str)
parser.add_argument('target_profile', action="store", type=str)
parser.add_argument('outfile', action="store", type=str)
parser.add_argument('join_on', action="store", type=str)


inputs = parser.parse_args()

make_labels(inputs.annotation_file, inputs.target_profile, inputs.outfile, inputs.join_on)






    
