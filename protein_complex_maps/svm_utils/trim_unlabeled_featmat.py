from __future__ import print_function
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it


def main():

    parser = argparse.ArgumentParser(description="Remove unlabeled entries from feature matrix")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                    help="Filename of output")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else 1st column")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)

    #kdrew: trim feature table to just the label column and the specified feature columns
    feature_table_trim = feature_table.query("%s != 0" % args.label_column)

    feature_table_trim.to_csv(args.output_filename, sep=args.sep)


if __name__ == "__main__":
    main()


