
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
    parser.add_argument("--global_preset", dest="global_preset", action = "store_true", default = False,required = False, 
                                    help = "If choosing custom ranges to normalize between, default False")
    parser.add_argument("--max", dest="maximum", type = float, action = "store", required = False, 
                                    help = "Choose preset maximum")
    parser.add_argument("--min", dest="minimum", type = float, action = "store", required = False, 
                                    help = "Choose preset minimum value")

    parser.add_argument("--inverse", action="store_true", dest="inverse", required=False, default=False,
                                    help="If True, rescale 0-n to 1-0, so that minimum value of feature is zero, ex. for euclidean distances, default=True")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)
   
    #Format features 
    features = read_feature_sel(args.features, args.feature_file, feature_table)

    if args.global_preset:
        if args.maximum:
            preset_maximum = args.maximum

        else:
            preset_maximum = feature_table[features].max().values[0]
            print(preset_maximum)


 
        if args.minimum:
             preset_minimum = args.minimum

        else:
             preset_minimum = feature_table[features].min().values[0]
             print(preset_minimum)

        feature_table = normalize_features_preset(features, feature_table, preset_maximum, preset_minimum, args.inverse)
    #Normalize selected features
    else:
        feature_table = normalize_features_col(features, feature_table, args.inverse)


    feature_table.to_csv(args.output_filename, index=False)


def normalize_features_preset(features, feature_table, preset_maximum, preset_minimum, inverse=True):
    '''
    Takes either an nargs+ list of features, or a file of feature names, one per line
    '''

    #Normalize each column to fixed values 
    if inverse == True:
        feature_table_norm = feature_table[features].apply(lambda x: 1- (x - preset_minimum)/(preset_maximum - preset_minimum), axis=0)

    else:
        feature_table_norm = feature_table[features].apply(lambda x: (x - preset_minimum)/(preset_maximum - preset_minimum), axis=0)

    feature_table[features] = feature_table_norm

    return feature_table
    

def normalize_features_col(features, feature_table, inverse=True):
    #Calculate column wise
    if inverse == True:
        feature_table_norm = feature_table[features].apply(lambda x: 1-  (x - x.min()) / (x.max() - x.min()), axis = 0)

    else:
        feature_table_norm = feature_table[features].apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)

    feature_table[features] = feature_table_norm

    return feature_table
 


def read_feature_sel(features, feature_file, feature_table):
    '''
    Takes either an nargs+ list of features, or a file of feature names, one per line
    Formats into list of features
    If none provided, use the features from the feature matrix
    '''
    if feature_file:
         with open(feature_file, "r") as featfile: 
             features = [x for x in featfile.read().split('\n') if x != '']

    elif features:
         features = features
    
    else: 
        #If no feature selection, keep all columns
        print("Warning: using all columns as features") 
        features  = feature_table.columns

    print("Normalizing the following features: ", features)
    return features


if __name__ == "__main__":
    main()
