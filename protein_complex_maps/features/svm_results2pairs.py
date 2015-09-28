
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu


def main():

    parser = argparse.ArgumentParser(description="Merges results from svm with ids of ppi pairs")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix, (full labeled matrix, uses only 0-labeled)")
    parser.add_argument("--input_results", action="store", dest="input_results", required=True, 
                                    help="Filename of results from libsvm with probability")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading feature csv, default=$")
    parser.add_argument("--add_prob", action="store_true", dest="add_prob", required=False, default=False,
                                    help="Flag for adding SVM probability of being true, default=False")
    parser.add_argument("--label_not0", action="store_true", dest="label_not0", required=False, default=False,
                                    help="Flag for only outputting rows with labels other than 0 (usually -1, 1), default=False")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")

    args = parser.parse_args()

    results = pd.read_csv(args.input_results,sep=' ')
    results.columns = ['svm_predicted_label','svm_pos_prob','svm_neg_prob']

    features = pd.read_csv(args.feature_matrix,sep=args.sep)
    if args.label_not0:
        only_features = features[features['label'] != 0]
    else:
        only_features = features[features['label'] == 0]
    only_features.index = range(len(only_features))

    only_features_results = pd.concat([results,only_features],axis=1)

    only_features_results_sorted = only_features_results.sort('svm_pos_prob',ascending=False)
    only_features_results_sorted_key2int = only_features_results_sorted[['key1','key2']].apply(lambda row: map(int,row))
    only_features_results_sorted_key2int_noself = only_features_results_sorted_key2int[only_features_results_sorted_key2int['key1'] != only_features_results_sorted_key2int['key2']]

    if args.add_prob:
        only_features_results_sorted_key2int.columns = ['key1_int','key2_int'] 
        #kdrew: adding full feature and results sorted (to get a handle on probability) 
        features_results_sorted_key2int = pd.concat([only_features_results_sorted_key2int,only_features_results_sorted], axis=1)
        #kdrew: remove self
        features_results_sorted_key2int_noself = features_results_sorted_key2int[features_results_sorted_key2int['key1_int'] != features_results_sorted_key2int['key2_int']]
        #kdrew: trim to just keys and prob
        keys_wprob = features_results_sorted_key2int_noself[['key1_int','key2_int','svm_pos_prob']]
        #kdrew: remove duplicate rows where the keys are the same
        keys_nodups = keys_wprob.drop_duplicates(['key1_int','key2_int'])
        print keys_nodups
    else:
        keys_only = only_features_results_sorted_key2int_noself[['key1','key2']]
        keys_nodups = keys_only.drop_duplicates()

    keys_nodups.to_csv(args.output_file,sep="\t",index=False,header=False)



if __name__ == "__main__":
    main()


