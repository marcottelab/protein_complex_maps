from __future__ import print_function
import pandas as pd
import argparse

# This is an all purpose script for annotating table with new columns
# 
#

def make_labels(annotation, infile, outfile, join_on, sep1, sep2, suffix, outfile_sep, target_col):

   #output from get_elution_ids.py
  
   annot=pd.read_csv(annotation, index_col=False, sep=sep1)
   print(annot) 

   annot = annot.set_index([join_on])

   if target_col:
      annot = annot[[target_col]]


   prof = pd.read_csv(infile, index_col=False, sep=sep2)
   print(prof)
   prof = prof.set_index([join_on])


   annotated = prof.join(annot, how="left", rsuffix=suffix)
   #print annotated
   annotated.to_csv(outfile, sep=str(outfile_sep))

   
parser = argparse.ArgumentParser(description='Short sample app')

parser.add_argument('annotation_file', action="store", type=str, help="table with annotations")
parser.add_argument('infile', action="store", type=str, help="Table to be annotated")
parser.add_argument('outfile', action="store", type=str)
parser.add_argument('join_on', action="store", type=str, help="common colname to join two files")
parser.add_argument('sep1', action="store", type=str, help="sep for annot file")
parser.add_argument('sep2', action="store", type=str, help="sep for target")
parser.add_argument('--suffix', action="store", type=str, default = "_annot", required=False, help="If annotation and infiles have common header names, requires suffix to join")
parser.add_argument('--outfile_sep', action="store", type=str, default=",", required=False, help="Separator for the created annotated file")
parser.add_argument('--target_col', action="store", dest="target_col",  required=False, default=None, help="designate a specific column to full from annotation file")


inputs = parser.parse_args()

make_labels(inputs.annotation_file, inputs.infile, inputs.outfile, inputs.join_on, inputs.sep1, inputs.sep2, inputs.suffix, inputs.outfile_sep, inputs.target_col)






    
