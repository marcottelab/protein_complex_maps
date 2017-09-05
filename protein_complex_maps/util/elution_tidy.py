
import sys
import numpy as np
import argparse
import pickle
import pandas as pd
import gc
from scipy import stats
#from guppy import hpy

def main():

    parser = argparse.ArgumentParser(description="Converts elution profile into tidy format")
    parser.add_argument("--input_elutions", action="store", dest="input_elutions", nargs='+', required=True, 
                                            help="Filename of elution profiles")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True,
                                            help="Filename of output file")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=',', 
                                            help="Separator to use when parsing elution profile")
    parser.add_argument("--calc_percentile", action="store_true", dest="calc_percentile", required=False, default=False, 
                                            help="Calculate percentiles")
    args = parser.parse_args()


    combined_melt_sparse = pd.DataFrame()
    for elution_fname in args.input_elutions:
        df = pd.read_csv(elution_fname,sep=args.sep)
        df.columns = ['prot_id'] + [x for x in df.columns[1:]]
        #print df.index
        #print df
        #kdrew: TODO check to see if TotalCount is in the header and if so ignore that column
        ##kdrew: makes full matrix = tidy, match ProtID with fraction names, start at column 2 because 0 and 1 are ProtID and TotalCounts
        df_melt = pd.melt(df, id_vars=df.columns[0], value_vars=[x for x in df.columns[1:]], var_name="fraction",value_name="abundance")
        df_melt_sparse = df_melt[df_melt.abundance != 0.0]
        combined_melt_sparse = pd.concat([combined_melt_sparse,df_melt_sparse])
        if args.calc_percentile:
            combined_melt_sparse['abundance_percentile'] = combined_melt_sparse['abundance'].apply(lambda x: stats.percentileofscore(combined_melt_sparse['abundance'].values,x))
        #print df_melt_sparse


    combined_melt_sparse.to_csv(args.out_filename)
    print combined_melt_sparse

    



if __name__ == "__main__":
    main()


