import pandas as pd
import sys
import argparse


def make_wide(identified_elution, orthology_file, protein_wideform, raw_loc):

   #output from get_elution_ids.py
   elut=pd.read_csv(identified_elution, index_col=False)
   
   elut= elut[['FractionID', 'GroupID', 'Total_SpecCounts']]

   #changing from long to wide format for elution profiles
   wide = elut.pivot(index='GroupID', columns = 'FractionID', values='Total_SpecCounts')  
   #print wide

   wide = wide.fillna(0)


   raw_outfile = raw_loc + identified_elution.replace("_elution_", "_raw_wide_elution_")
   wide.to_csv(raw_outfile)

   wide['Total'] = wide.sum(axis=1)


   #saving and importing csv fixes header columns from pivot format
   wide2 = pd.read_csv(raw_outfile)
   wide2 = wide2.set_index(['GroupID'])
 

   #Pull spectral count mapped to individual proteins
   protein_wide = pd.read_csv(protein_wideform, index_col=False)
   protein_wide = protein_wide[['ProteinID', 'Total']]
   protein_wide['Total'] = protein_wide['Total'].astype(str)

   #print protein_wide.loc[protein_wide['ProteinID'] == "tr|W5DIX6|W5DIX6_WHEAT"]

   protein_wide['Total'] = protein_wide['Total'].str.replace("\.0", "") #To remove decimal, everything is an integer
   #print protein_wide.loc[protein_wide['ProteinID'] == "tr|W5DIX6|W5DIX6_WHEAT"]

   protein_wide['ProteinTotals'] = protein_wide[['ProteinID', 'Total']].apply(lambda x:":".join(x), axis=1)   
   protein_wide = protein_wide[['ProteinID', 'ProteinTotals']]
   #print protein_wide.loc[protein_wide['ProteinID'] == "tr|W5DIX6|W5DIX6_WHEAT"]
   protein_wide = protein_wide.set_index(['ProteinID'])
#   print protein_wide


   #Make in to table that is ProteinID ProteinID:Total
   #Join to annot

   #Pull annotations from orthology file (eggnog_output)
   annot = pd.read_csv(orthology_file, sep="\t")
   annot = annot[['GroupID', 'ProteinID', 'Annotation']]
   #annot = annot.set_index(['GroupID']) 
 
   annot2 = annot[['GroupID', 'ProteinID']]
   annot2 = annot2.set_index(['ProteinID'])

#   print annot2
#   print protein_wide
   annot3 = annot2.join(protein_wide, how="left")
   annot3 = annot3.reset_index()
  
   annot3 = annot3[['GroupID', 'ProteinTotals']]
#   print annot3
   abundance =   annot3.groupby(['GroupID'])['ProteinTotals'].apply(lambda x: ' '.join(map(str, x))).reset_index()
   abundance['ProteinTotals'] = abundance['ProteinTotals'].str.replace("nan ", "")
   abundance['ProteinTotals'] = abundance['ProteinTotals'].str.replace(" nan", "")
   abundance = abundance.set_index(['GroupID'])


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
   #print alt_wide
   #print alt_wide

   alt_wide_abundance = alt_wide.join(abundance, how = "left")

   alt_outfile = identified_elution.replace("_elution_", "_alt_wide_elution_")
   print alt_wide_abundance
   alt_wide_abundance.to_csv(alt_outfile)

    
parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('identified_elution', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
parser.add_argument('protein_wide', action="store", type=str)
parser.add_argument('raw_outfile_loc', action="store", type=str)
inputs = parser.parse_args()

make_wide(inputs.identified_elution, inputs.orthology_file, inputs.protein_wide, inputs.raw_outfile_loc)






    
