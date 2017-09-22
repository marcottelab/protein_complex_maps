from __future__ import print_function
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it


def main():

    parser = argparse.ArgumentParser(description="Converts a dataframe of features to libsvm format")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                    help="Filename of output")
    parser.add_argument("--features", action="store", dest="features", nargs='+', required=False, 
                                    help="Names of features to output")
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False, 
                                    help="File containing features, one per line")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else 1st column")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)
    if args.feature_file:
         with open(args.feature_file, "r") as featfile: 
             features = [x for x in featfile.read().split('\n') if x != '']

         print(features)

    elif args.features:
         features = args.features
    
    else: 
        #If no feature selection, keep all columns
        print("Warning: using all columns as features") 
        features  = feature_table.columns

    #kdrew: trim feature table to just the label column and the specified feature columns
    feature_table_trim = feature_table[[args.label_column] + features]


    outfile = open(args.output_filename,"wb")
    for i in feature_table_trim.to_records(index=False):
        #kdrew: build libsvm formatted string, [label 1 or -1 or 0] 1:[feature value] 2:[feature value] ...

        #cdm for large feature mats, write as you go instead concatenating
        output_string = "%s %s\n" % (int(i[0]), " ".join(["%s:%s" % (j+1,x) for j, x in enumerate(i.tolist()[1:])]))
        outfile.write(output_string)


    outfile.close()



if __name__ == "__main__":
    main()


