from __future__ import print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

def main():

    parser = argparse.ArgumentParser(description="Selects multiiple selections of columns from a list")
    parser.add_argument("--libsvm1_scaled", action="store", dest= "df1", required=True)
    parser.add_argument("--feature_columns", action="store", required=True)
    parser.add_argument("--feature_header", action="store", required=True)
    args = parser.parse_args()
 
    keep = ['value']

    with open(args.feature_header, "r") as fh:
        feature_tmp = fh.read().rstrip()
        print(feature_tmp)
        feature_names = feature_tmp.split(",")
    head= keep + feature_names
    print(head)
 
    print("opening df1")
    df1 = pd.read_table(args.df1, sep=" ", header=None)
    print(df1)
    df1.drop(df1.columns[len(df1.columns)-1], axis=1, inplace=True)
    print(df1)
    df1.columns = head



    print(df1)

    with open(args.feature_columns, "r") as fc:
         f_files = [line.rstrip() for line in fc]
         for f_file in f_files:
             with open(f_file, "r") as ff:
                 features = [line.rstrip() for line in ff]        
                 features = keep + features
                 print(df1[features])
                 outdf = df1[features]
                 outfile = f_file.replace(".csv", "") + "_" + args.df1     
                 outdf.to_csv(outfile, index=False, header=False, sep=" ")
 

#    print(df1)
#    df1['pairset'] = map(frozenset, zip( df1['ID1'].values, df1['ID2'].values ))
#    print("df1 index setting")
#    df1 = df1.set_index(['pairset'])
#    print("opening df2")
#    df2 = pd.read_table(args.df2, sep=",")
#    print(df2)
#    df2['pairset'] = map(frozenset, zip( df2['ID1'].values, df2['ID2'].values ))
#    print("df2 index set")
#
#    df1 = df1.drop(['ID1', 'ID2'], 1)
#    df2 = df2.drop(['ID1', 'ID2'], 1)
#
#    dfall = df1.join(df2, how='outer')
#
#

    #dfall = dfall.reset_index()
    #dfall['ID1'] = [list(x)[0] for x in dfall['pairset'].values]
    #dfall['ID2'] = [list(x)[1] for x in dfall['pairset'].values]

 
    #dfall = dfall.drop('pairset', 1)

    #dfall.to_csv( args.out_filename, index=False)
    



if __name__ == "__main__":
    main()


