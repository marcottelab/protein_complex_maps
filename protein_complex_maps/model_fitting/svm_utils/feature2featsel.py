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
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False,
                                    help="File containing features, one per line")

    parser.add_argument("--features", action="store", dest="features", nargs='+', required=False, 
                                    help="Names of features to output")
    parser.add_argument("--id_column", action="store", dest="id_column", required=False, default='ID', 
                                    help="Name of ID column, default='ID' if present, else 1st column")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else last column")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)
    #kdrew: trim feature table to just the label column and the specified feature columns

    if args.feature_file:
         with open(args.feature_file, "r") as featfile:
             features = [x for x in featfile.read().split('\n') if x != '']

    elif args.features:
         features = args.features
    else:
         print("Need to provide features to select")
         return(0)
    print(features)

    feature_table_trim = feature_table[[args.id_column] + features + [args.label_column]]

    
    feature_table_trim.to_csv(args.output_filename, index=False)


if __name__ == "__main__":
    main()


