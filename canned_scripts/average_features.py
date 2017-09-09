#! /usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Average a list of feature files")
parser.add_argument("infiles", nargs="+",help='Feature files. Optionally takes a wildcard like *.csv')
parser.add_argument("outfile",help="Name of file to write to")
parser.add_argument("--average_type", default="mean", choices=["mean","median"], help="Choose mean or median values")
parser.add_argument("--retain_columns",action='store_true',help="Retain the input columns in the ouput")
parser.add_argument("--log",action='store_true',help="Write log of inputs")
args = parser.parse_args()

if __name__ == '__main__':

    files_processed = []
    is_first = True
    for ind,f in enumerate(args.infiles):
        print "Reading in {}".format(f)
        df = pd.read_csv(f)
        feature = df.columns[2]
        df.columns = ["ID1","ID2","file"+str(ind)]
        if is_first:
            merged = df
            is_first = False
            continue
        merged = merged.merge(df,on=["ID1","ID2"],how='outer')
        files_processed.append(f)

    if args.average_type == "mean":
        merged["mean_{}".format(feature)] = merged.iloc[:,2:].mean(axis=1)
    else:
        merged["median_{}".format(feature)] = merged.iloc[:,2:].median(axis=1)
        
    if args.retain_columns == False:
        drop_cols = [col for col in merged.columns if "file" in col]
        merged.drop(drop_cols,axis=1,inplace=True)
    merged.to_csv(args.outfile,index=False)

    if args.log:
        with open("average_features.log",'w') as out:
            _ = [out.write("Processed {}".format(i) + "\n") for i in files_processed]

