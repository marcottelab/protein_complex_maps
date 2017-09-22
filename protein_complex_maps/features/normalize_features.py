from __future__ import print_function
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it


def main():

    parser = argparse.ArgumentParser(description="Normalize selected columns of a feature matrix")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                    help="Filename of output")
    parser.add_argument("--features", action="store", dest="features", nargs='+', required=False, 
                                    help="Names of features to output")
    parser.add_argument("--feature_file", action="store", dest="feature_file", required=False, 
                                    help="File containing features, one per line")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")
    parser.add_argument("--inverse", action="store_true", dest="inverse", required=False, default=True,
                                    help="If True, rescale 0-n to 1-0, so that minimum value of feature is zero, ex. for euclidean distances, default=True")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)
   
    #Format features 
    features = read_feature_sel(args.features, args.feature_file, feature_table)

    #Normalize selected features
    feature_table = normalize_features(features, feature_table, args.inverse)

    feature_table.to_csv(args.output_filename, index=False)


def normalize_features(features, feature_table, inverse=True):
    '''
    Takes either an nargs+ list of features, or a file of feature names, one per line
    '''

    #cdm: trim feature table to just the specified feature columns
    feature_table_trim = feature_table[features]

    if inverse == True:
        feature_table_norm = feature_table_trim.apply(lambda x: 1- (x/x.max()), axis=0)

    else:
        feature_table_norm = feature_table_trim.apply(lambda x: (x/x.max()), axis=0)

    feature_table[features] = feature_table_norm

    return feature_table



def read_feature_sel(features, feature_file, feature_table):
    '''
    Takes either an nargs+ list of features, or a file of feature names, one per line
    Formats into list of features
    If none provided, use the features from the feature matrix
    '''
    if feature_file:
         with open(args.feature_file, "r") as featfile: 
             features = [x for x in featfile.read().split('\n') if x != '']

    elif features:
         features = args.features
    
    else: 
        #If no feature selection, keep all columns
        print("Warning: using all columns as features") 
        features  = feature_table.columns

    print("Normalizing the following features: ", features)
    return features


if __name__ == "__main__":
    main()


