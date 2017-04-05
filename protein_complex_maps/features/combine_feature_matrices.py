from __future__ import print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

#This needs to be updated to not use frozenset

def main():

    parser = argparse.ArgumentParser(description="Builds feature matrix by concatenating pairs files")
    parser.add_argument("--next_feature_matrix", action="store", dest= "df1", required=True)
    parser.add_argument("--prev_feature_matrix", action="store", dest="df2", required=False)
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True, default=None, 
                                    help="Filename of output file")
    args = parser.parse_args()
    
  
    print("opening df1")
    df1 = pd.read_table(args.df1, sep=",")
    print(df1)
    df1['pairset'] = map(frozenset, zip( df1['ID1'].values, df1['ID2'].values ))
    print("df1 index setting")
    df1 = df1.set_index(['pairset'])
    print("opening df2")
    df2 = pd.read_table(args.df2, sep=",")
    print(df2)
    df2['pairset'] = map(frozenset, zip( df2['ID1'].values, df2['ID2'].values ))
    print("df2 index set")

    df1 = df1.drop(['ID1', 'ID2'], 1)
    df2 = df2.drop(['ID1', 'ID2'], 1)

    dfall = df1.join(df2, how='outer')



    dfall = dfall.reset_index()
    dfall['ID1'] = [list(x)[0] for x in dfall['pairset'].values]
    dfall['ID2'] = [list(x)[1] for x in dfall['pairset'].values]

 
    dfall = dfall.drop('pairset', 1)

    dfall.to_csv( args.out_filename, index=False)
    



if __name__ == "__main__":
    main()


