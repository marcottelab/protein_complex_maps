from __future__ import print_function
import pandas as pd
import argparse

# This is an all purpose script for annotating table with new columns
# 
#

def do_annotation(annotation, infile, outfile, join_on, jointype, sep1, sep2, suffix, outfile_sep, target_col):

   #output from get_elution_ids.py
  
   annot=pd.read_csv(annotation, index_col=False, sep=sep1)
   print(annot) 

   print(join_on)
   annot = annot.set_index(join_on)
   if target_col:
      annot = annot[[target_col]]


   tbl = pd.read_csv(infile, index_col=False, sep=sep2)
   print(tbl)
   print(join_on)
   tbl = tbl.set_index(join_on)


   annotated = tbl.join(annot, how=jointype, rsuffix=suffix)
   annotated.to_csv(outfile, sep=str(outfile_sep))

   
parser = argparse.ArgumentParser(description='Very general script for annotating a table with column[s] from a different table')

parser.add_argument('-a', '--annotation_file', action="store", type=str, required=True,  help="table with annotations")
parser.add_argument('-i','--infile', action="store", type=str, required=True, help="Table to be annotated")
parser.add_argument('-o','--outfile', action="store", type=str, required=True)
parser.add_argument('-jo','--join_on', action="store", type=str, nargs='+', required=True, help="common colname to join two files")
parser.add_argument('-jt','--join_type', action="store", type=str, required=False, default="left", choices=["left", "right", "inner", "outer"], help="Type of join. [left, right, inner, outer]")
parser.add_argument('-s1','--sep1', action="store", type=str, default=',', help="sep for annot file")
parser.add_argument('-s2','--sep2', action="store", type=str, default=',', help="sep for target")
parser.add_argument('-sf','--suffix', action="store", type=str, default = "_annot", required=False, help="If annotation and infiles have common header names, requires suffix to join")
parser.add_argument('-os','--outfile_sep', action="store", type=str, default=",", required=False, help="Separator for the created annotated file")
parser.add_argument('-t','--target_col', action="store", dest="target_col",  required=False, default=None, help="designate a specific column to full from annotation file")


inputs = parser.parse_args()

do_annotation(inputs.annotation_file, inputs.infile, inputs.outfile, inputs.join_on, inputs.join_type, inputs.sep1, inputs.sep2, inputs.suffix, inputs.outfile_sep, inputs.target_col)
