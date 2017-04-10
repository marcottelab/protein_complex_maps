import pandas as pd
import sys
import argparse
import logging

def make_wide(identified_elution, orthology_file):


   #output from get_elution_ids.py
   elut=pd.read_csv(identified_elution, index_col=False)
   
   elut= elut[['FractionID', 'GroupID', 'Total_SpecCounts']]

   #changing from long to wide format for elution profiles
   wide = elut.pivot(index='GroupID', columns = 'FractionID', values='Total_SpecCounts')  
   #print wide

   wide = wide.fillna(0)


   raw_outfile = identified_elution.replace("_elution_", "_raw_wide_elution_")
   wide.to_csv(raw_outfile)

   wide['Total'] = wide.sum(axis=1)


   #saving and importing csv fixes header columns from pivot format
   wide2 = pd.read_csv(raw_outfile)
   wide2 = wide2.set_index(['GroupID'])

   #Pull annotations from orthology file (eggnog_output)
   annot = pd.read_csv(orthology_file, sep="\t")
   annot = annot[['GroupID', 'ProteinID', 'Annotation']]

   #put one protein per line
   annotmulti = annot.set_index(['GroupID']) 

   multirow_wide = annotmulti.join(wide, how = "inner")

   multirow_outfile = identified_elution.replace("_elution_", "_multirow_wide_elution_")
    
   multirow_wide.to_csv(multirow_outfile)

   #Consolidate proteins from a group into one row
   alt_wide_labels = annot.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()

   #Get annotations and 
   #ungrouped_alt_wide = annot.join(wide, how = "right")

   #ungrouped_alt_wide = ungrouped_alt_wide.reset_index()
   #http://stackoverflow.com/questions/27298178/concatenate-strings-from-several-rows-using-pandas-groupby
   #alt_wide_labels = ungrouped_alt_wide.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
#   print grouped_elut

   alt_wide_labels = alt_wide_labels.set_index(["GroupID"])
   alt_wide = alt_wide_labels.join(wide, how = "inner")
   print alt_wide
   #print alt_wide

   alt_outfile = identified_elution.replace("_elution_", "_alt_wide_elution_")
   alt_wide.to_csv(alt_outfile)

    
parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('identified_elution', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
inputs = parser.parse_args()

make_wide(inputs.identified_elution, inputs.orthology_file)






    
