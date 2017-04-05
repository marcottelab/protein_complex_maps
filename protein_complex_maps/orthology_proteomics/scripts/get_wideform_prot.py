import pandas as pd
import sys
import argparse


def make_wide(identified_elution):

   #output from get_elution_ids.py
   elut=pd.read_csv(identified_elution, index_col=False)
   
   elut= elut[['FractionID', 'ID', 'Total_SpecCounts']]

   #changing from long to wide format for elution profiles
   wide = elut.pivot(index='ID', columns = 'FractionID', values='Total_SpecCounts')  
   #print wide

   wide = wide.fillna(0)


   raw_outfile = identified_elution.replace("_elution_", "_raw_wide_elution_")
   wide.to_csv(raw_outfile)

   wide['Total'] = wide.sum(axis=1)


   alt_outfile = identified_elution.replace("_elution_", "_alt_wide_elution_")
   wide.to_csv(alt_outfile)

    
parser = argparse.ArgumentParser(description='Get wide form')

parser.add_argument('identified_elution', action="store", type=str)
inputs = parser.parse_args()

make_wide(inputs.identified_elution)






    
