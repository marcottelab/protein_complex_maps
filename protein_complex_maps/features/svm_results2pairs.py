
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
                                    help="Filename of input feature matrix")
    parser.add_argument("--input_results", action="store", dest="input_results", required=True, 
                                    help="Filename of results from libsvm with probability")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading feature csv, default=$")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")

    args = parser.parse_args()

    #only0_results = pd.read_csv('/home/kdrew/data/protein_complex_maps/v35_features/blake_bioplex_feature_merge_wkeys_labeled_nan2num.libsvm0.scale.resultsWprob',sep=' ')
    only0_results = pd.read_csv(args.input_results,sep=' ')
    only0_results.columns = ['svm_predicted_label','svm_pos_prob','svm_neg_prob']

    #bf_bioplex_feature_merge_labeled_nan2num.to_csv("/home/kdrew/data/protein_complex_maps/v35_features/blake_bioplex_feature_merge_wkeys_labeled_nan2num.txt",sep="$")
    features = pd.read_csv(args.feature_matrix,sep=args.sep)
    only0_features = features[features['label'] == 0]
    only0_features.index = range(len(only0_features))

    only0_features_results = pd.concat([only0_results,only0_features],axis=1)

    only0_features_results_sorted = only0_features_results.sort('svm_pos_prob',ascending=False)
    only0_features_results_sorted_key2int = only0_features_results_sorted[['key1','key2']].apply(lambda row: map(int,row))
    only0_features_results_sorted_key2int_noself = only0_features_results_sorted_key2int[only0_features_results_sorted_key2int['key1'] != only0_features_results_sorted_key2int['key2']]
    #only0_features_results_sorted_key2int_noself[['key1','key2']].to_csv("/home/kdrew/data/protein_complex_maps/v35_features/blake_bioplex_feature_merge_wkeys_labeled_nan2num.libsvm0.scale.features_resultsWprob_pairs_noself.txt",sep="\t",index=False,header=False)
    keys_only = only0_features_results_sorted_key2int_noself[['key1','key2']]
    keys_nodups = keys_only.drop_duplicates()
    #only0_features_results_sorted_key2int_noself[['key1','key2']].to_csv(args.output_file,sep="\t",index=False,header=False)
    keys_nodups.to_csv(args.output_file,sep="\t",index=False,header=False)



if __name__ == "__main__":
    main()


