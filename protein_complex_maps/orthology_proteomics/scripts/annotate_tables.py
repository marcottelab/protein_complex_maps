from __future__ import print_function
import pandas as pd
import sys
import argparse


def make_labels(annotation, profile, outfile, join_on, sep1, sep2, suffix, outfile_sep, target_col):

   #output from get_elution_ids.py
  
   annot=pd.read_csv(annotation, index_col=False, sep=sep1)
   #print annot 

   annot = annot.set_index([join_on])

   if target_col:
      annot = annot[[target_col]]


   prof = pd.read_csv(profile, index_col=False, sep=sep2)
   #print prof
   prof = prof.set_index([join_on])


   annotated = prof.join(annot, how="left", rsuffix=suffix)
   #print annotated
   annotated.to_csv(outfile, sep=str(outfile_sep))

   
parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('annotation_file', action="store", type=str)
parser.add_argument('target_profile', action="store", type=str)
parser.add_argument('outfile', action="store", type=str)
parser.add_argument('join_on', action="store", type=str)
parser.add_argument('sep1', action="store", type=str)
parser.add_argument('sep2', action="store", type=str)
parser.add_argument('suffix', action="store", type=str, default = "")
parser.add_argument('outfile_sep', action="store", type=str)
parser.add_argument('--target_col', action="store", dest="target_col",  required=False, default=None)


inputs = parser.parse_args()

make_labels(inputs.annotation_file, inputs.target_profile, inputs.outfile, inputs.join_on, inputs.sep1, inputs.sep2, inputs.suffix, inputs.outfile_sep, inputs.target_col)






    
