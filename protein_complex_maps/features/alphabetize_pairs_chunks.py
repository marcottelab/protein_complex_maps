from __future__ import print_function
import os
import argparse
import pandas as pd
import numpy as np
import timeit

def alphabetized_check(df, column_ids, sample_size=1000):
    # Fails with Pandas 1.0+
    if sample_size > len(df):
       sample_size = int(len(df)/10)
    print(column_ids)
    #kdrew: sample a portion of the passed in dataframe
    df_sample = df.sample(sample_size)
    #kdrew: alphabetize the passed in columns and record the first column in 'alpha_id1'
    df_sample['alpha_id1'] = df_sample[column_ids].apply(sorted,axis=1)[column_ids[0]]
    #kdrew: check to see if all of the first alphabetized ids are the same as the first passed in id column
    ret_value = all(df_sample['alpha_id1'] == df_sample[column_ids[0]])
    return ret_value


def main():
    '''
    This function takes a file of pairwise ID features:
    A B value
    D C value
    And orders the identifiers:
    A B value
    C D value
    This is because different feature makers may order the pairwise IDs differently
    '''


    parser = argparse.ArgumentParser(description="Selects multiple selections of columns from a list")
    parser.add_argument("--feature_pairs", action="store", dest= "df", required=True)
    parser.add_argument("--outfile", action="store", required=True, 
                                    help="Filename of comma separated output")
    parser.add_argument("--sep", action = "store",dest='sep' , default=' ', required=False)
    parser.add_argument("--columns", nargs='+', action="store", type=int, dest="columns2alphabetize", required=False, default=[0,1], 
                                    help="List of columns (positions) to alphabetize, default = 0,1")
    parser.add_argument("--header", action="store", dest="header", required=False, default=True, 
                                    help="Do input pairs files have headers?")
    parser.add_argument("--chunksize", action="store", dest="chunksize", type=int, required=False, default=100000, 
                                    help="Chunksize for processing, default = 100000")


    args = parser.parse_args()

    if args.sep == '\\t':
           args.sep = '\t'
 
    if args.header == True:
        iterator = pd.read_csv(args.df, sep=args.sep, header=0, engine="python",  chunksize=args.chunksize, iterator = True)
    if args.header == False:
        iterator  = pd.read_csv(args.df, sep=args.sep, header=None, engine="python", chunksize=args.chunksize, iterator = True)


    firstpass = True
    for df in iterator:
        alphabetized_df = alphabetize_df(df, args.columns2alphabetize)

        if firstpass == True:
            if args.header == True:
                alphabetized_df.to_csv(args.outfile, header=True, index=False, sep=args.sep)
            if args.header == False:
                alphabetized_df.to_csv(args.outfile, header=False, index=False, sep=args.sep)

        else: # append without writing the header
            alphabetized_df.to_csv(args.outfile, mode='a', index=False, header=False, sep=args.sep)
        firstpass = False
  

def alphabetize_df(df, columns2alphabetize):

    try:
       intermediate_df =  df[df.columns[columns2alphabetize]].apply(sorted,axis=1, broadcast = True)

    except Exception as E:
       #Pandas 0.18 doesn't have the broadcast option, but returns dataframe by default from this command
       print(E)
       print("Exception in apply, trying alternate way for pre-broadcast pandas installations")
       intermediate_df =  df[df.columns[columns2alphabetize]].apply(sorted,axis=1)
   
    df[df.columns[columns2alphabetize]] = intermediate_df
   
    return df




if __name__ == "__main__":
    main()


