from __future__ import print_function
import pandas as pd
import sys
import argparse


def make_labels(kog_file, orthology_file, outfile):

   #output from get_elution_ids.py
   egg=pd.read_csv(orthology_file, index_col=False, sep="\t")

   kog=pd.read_csv(kog_file, index_col=False, sep="\t")
   print(egg)
   print(kog)  
   test1 = egg[(egg["GroupID"] == "KOG0002")]
   egg = egg.set_index(['GroupID'])
   kog = kog.set_index(['GroupID'])
 
  

   print(test1)


   egg.update(kog[['Annotation']])

   egg = egg.reset_index()
   test2 = egg[(egg["GroupID"] == "KOG0002")]
  

   print(test2)
   egg.to_csv(outfile, sep="\t", index=False)
    
parser = argparse.ArgumentParser(description='Update orthology annotations with new KOG standard')

#parser.add_argument('identified_eggion', action="store", type=str)
parser.add_argument('kog_file', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
#parser.add_argument('target_profile', action="store", type=str)
parser.add_argument('outfile', action="store", type=str)
inputs = parser.parse_args()

make_labels(inputs.kog_file, inputs.orthology_file, inputs.outfile)






    
