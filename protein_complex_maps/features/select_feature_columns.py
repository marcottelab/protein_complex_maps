from __future__ import print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

def main():

    parser = argparse.ArgumentParser(description="Selects multiple selections of columns from a list. This was for libsvm...")
    parser.add_argument("--libsvm1_scaled", action="store", dest= "df1", required=True)
    parser.add_argument("--feature_columns", action="store", required=True)
    parser.add_argument("--feature_header", action="store", required=True)
    args = parser.parse_args()
 
    keep = ['value']

    with open(args.feature_header, "r") as fh:
        feature_tmp = fh.read().rstrip()
        feature_names = feature_tmp.split(",")
    head= keep + feature_names
 
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
 


if __name__ == "__main__":
    main()


