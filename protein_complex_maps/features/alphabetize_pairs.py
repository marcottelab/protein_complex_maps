#from __future__ import print_function
import argparse
import pandas as pd
import numpy as np
import timeit
#from tqdm import tqdm
#python ../../scripts/alphabetize_pairs.py --feature_pairs arathtraesorysjbraolselml_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob_c32_g0078125_pairs_noself_nodups_wprob.txt --outfile old_arathtraesorysjbraolselml_interactions.tmp

def alphabetized_check(df, column_ids, sample_size=100):
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

    args = parser.parse_args()

    if args.header == True:
        df = pd.read_table(args.df, sep=args.sep, header=0)
    if args.header == False:
        df = pd.read_table(args.df, sep=args.sep, header=None)


    print(df)


    df = alphabetize_df(df, args.columns2alphabetize)
    if args.sep == '\\t':
           args.sep = '\t'
  
    if args.header == True:
        df.to_csv(args.outfile, header=True, index=False, sep=args.sep)
    if args.header == False:
        df.to_csv(args.outfile, header=False, index=False, sep=args.sep)


def alphabetize_df(df, columns2alphabetize):

    intermediate_df =  df[df.columns[columns2alphabetize]].apply(sorted,axis=1)
    print(intermediate_df)

    df[df.columns[columns2alphabetize]] = intermediate_df
   
    return df




if __name__ == "__main__":
    main()


