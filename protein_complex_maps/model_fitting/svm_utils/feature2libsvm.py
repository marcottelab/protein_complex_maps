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
                                    help="File containing one feature name per line")

    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else 1st column")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")

    args = parser.parse_args()

    if args.features:
         features = args.features

    elif args.feature_file:
        featfile = open(args.feature_file, "r")
        feats = featfile.read().split("\n")
        featfile.close()
        print(feats)
        features = ",".join(feats)
        print(features)

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)

    #kdrew: trim feature table to just the label column and the specified feature columns
    feature_table_trim = feature_table[[args.label_column] + args.features]

    output_string = ""
    for i in feature_table_trim.to_records(index=False):
        #kdrew: build libsvm formatted string, [label 1 or -1 or 0] 1:[feature value] 2:[feature value] ...
        output_string += "%s %s\n" % (int(i[0]), " ".join(["%s:%s" % (j+1,x) for j, x in enumerate(i.tolist()[1:])]))

    outfile = open(args.output_filename,"wb")
    outfile.write(output_string)
    outfile.close()



if __name__ == "__main__":
    main()


