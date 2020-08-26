#! /usr/bin/env python
from __future__ import print_function
import argparse
import functools
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Average a list of feature files")
parser.add_argument("infiles", nargs="+",help='Feature files. Optionally takes a wildcard like *.csv')
parser.add_argument("outfile",help="Name of file to write to")
parser.add_argument("--average_type", nargs="+", default="mean", choices=["mean","median",'max','min'], help="Choose mean or median values")
parser.add_argument("--retain_columns",action='store_true',help="Retain the input columns in the ouput")
parser.add_argument("--log",action='store_true',help="Write log of inputs")
parser.add_argument("--in_pickle",action='store_true',help="Whether to accept input as pickled DataFrame from pandas. Default `False`")
parser.add_argument("--out_pickle",action='store_true',help="Whether to save output as pickled DataFrame. Default `False`")
parser.add_argument("--z_transform",action='store_true',help="Perform Fisher's Z-tranform on correlation data before averaging")
parser.add_argument("--missing_penalization_factor",default=False,type=float,help="How strongly to penalize the absence of a pair from an elution. Larger = higher penalty")
args = parser.parse_args()

if __name__ == '__main__':

    files_processed = []
    is_first = True
    for ind,f in enumerate(args.infiles):
        print("Reading in {}".format(f))
        if args.in_pickle:
            df = pd.read_pickle(f)
            assert len(df.columns) == 3, "Pickled infiles should have 3 columns: ID1, ID2, feature"
            assert df.columns[0] == "ID1", "Pickled infiles should have 3 columns: ID1, ID2, feature"
            assert df.columns[1] == "ID2", "Pickled infiles should have 3 columns: ID1, ID2, feature"
            df.set_index(["ID1","ID2"],inplace=True)
        else:
            df = pd.read_csv(f,index_col=[0,1])
        assert len(df.columns) == 1, "Files should be DataFrames should have just one column"
        feature = df.columns[0]
        
        if "/" in f: # f contains a full path, probably
            df.columns = [f.split("/")[-1].split(".")[0]]
        else:
            df.columns = [f.split(".")[0]]
            
        if is_first:
            merged = df
            is_first = False
            continue
        merged = merged.join(df,how='outer')
        files_processed.append(f)
    
    if args.z_transform:
        merged[merged == 1.0] = .9999999999
        merged = np.arctanh(merged)
        
    ## To penalize missing data, add columns of zeros to the merged elutions
    ## The number of columns added is missing_penalization_factor * number_of_elutions.
    ## Idea is that adding a small number of columns with 0s will minimally affect the 
    ## that are widely observed and more strongly affect those found only in a few fractions
    ## So the best bet is to keep missing_penalization_factor low, maybe .2 or so.
    if args.missing_penalization_factor:
        n_to_add = int( args.missing_penalization_factor * len(merged.columns) )
        for i in xrange(n_to_add):
            merged["dummy{}".format(i)] = 0
        
    avgFuncD = {"mean": functools.partial( merged.mean ),
                "median": functools.partial( merged.median ),
                "max": functools.partial( merged.max ),
                "min": functools.partial( merged.min )}
    
    for avg in args.average_type:
        print("Performing {}".format(avg))
        merged["{}_{}".format(avg,feature)] = avgFuncD[avg](axis=1)
        
    if args.retain_columns == False:
        drop_cols = [col for col in merged.columns[:-len(args.average_type)]]
        merged.drop(drop_cols,axis=1,inplace=True)
    
    if args.out_pickle:
        merged.to_pickle(args.outfile)
    else:
        merged.to_csv(args.outfile)

    if args.log:
        with open("average_features.log",'w') as out:
            _ = [out.write("Processed {}".format(i) + "\n") for i in files_processed]

