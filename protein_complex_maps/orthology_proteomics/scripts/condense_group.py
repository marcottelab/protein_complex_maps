from __future__ import print_function
import pandas as pd
import argparse


def consolidate_rows(mapping_file, grouping_col, to_consolidate_col, outfile):

   
   map_df=pd.read_csv(mapping_file, index_col=False, sep="\t")
   
   
   final = map_df.groupby([grouping_col])[to_consolidate_col].apply(lambda x: ' '.join(x)).reset_index()

   final.to_csv(outfile, sep="\t", index_col=False)


parser = argparse.ArgumentParser(description='Take a two-column mapping and consolidate multiple entries into one row')
parser.add_argument('-m', 'mapping_file', action="store", type=str, help="a two column tsv")
parser.add_argument('-g', 'grouping_col', action="store", type=str, help="The column to group by")
parser.add_argument('-c', 'to_consolidate_col', action="store", type=str, help="the column to consolidate")
parser.add_argument('-o', 'outfile', action="store", type=str)
inputs = parser.parse_args()
consolidate_rows(inputs.mapping_file, inputs.grouping_col, inputs.to_consolidate_col, inputs.outfile)






    
