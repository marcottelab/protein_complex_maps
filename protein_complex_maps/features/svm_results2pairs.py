
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu

#kdrew: cannot have error checking in lambda function so create a function with error checking to be called by lambda
#kdrew: http://stackoverflow.com/questions/12451531/python-try-catch-block-inside-lambda
def tryconvertInt(value):
    #print "%s %s" % (value, type(value))
    ret = None
    try:
        ret = int(float(value))
        #print "no Error"
    except ValueError, TypeError:
        print "%s %s" % (value, type(value))
        print "Error"
        #return None
    return ret

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
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=False, default=['key1','key2'],
                                    help="List of columns that specify ids in feature matrix, default: key1 key2")
    parser.add_argument("--int_convert", action="store_true", dest="int_convert", required=False, default=False, 
                                    help="Convert id_columns to int, default=False")

    args = parser.parse_args()

    results = pd.read_csv(args.input_results,sep=' ')
    results.columns = ['svm_predicted_label','svm_pos_prob','svm_neg_prob']

    print results
    print " "

    features = pd.read_csv(args.feature_matrix,sep=args.sep)
    print features
    print " "
    if args.label_not0:
        only_features = features[features['label'] != 0]
    else:
        only_features = features[features['label'] == 0]
    only_features.index = range(len(only_features))
    print only_features
    print " "

    #kdrew: put the results and the features together
    only_features_results = pd.concat([results,only_features],axis=1)
    print only_features_results
    print " "

    #kdrew: cleaning up space
    del results
    del only_features

    only_features_results = only_features_results.astype('object')
    print only_features_results
    print " "

    only_features_results_sorted = only_features_results.sort('svm_pos_prob',ascending=False)
    print only_features_results_sorted
    print " "
    del only_features_results

    if args.int_convert:
        #only_features_results_sorted_key2int = only_features_results_sorted[['key1','key2']].apply(lambda row: map(tryconvertInt,row))
        only_features_results_sorted_keys = only_features_results_sorted[args.id_columns].apply(lambda row: map(tryconvertInt,row))
    else:
        only_features_results_sorted_keys = only_features_results_sorted[args.id_columns]

    print only_features_results_sorted_keys
    print " "

    #only_features_results_sorted_key2int = only_features_results_sorted_key2int.dropna().astype(object)
    only_features_results_sorted_keys_noself = only_features_results_sorted_keys[only_features_results_sorted_keys[args.id_columns[0]] != only_features_results_sorted_keys[args.id_columns[1]]]

    if args.add_prob:
        only_features_results_sorted.columns = ['key1_int','key2_int'] 
        #kdrew: adding full feature and results sorted (to get a handle on probability) 
        features_results_sorted_key2int = pd.concat([only_features_results_sorted_key2int,only_features_results_sorted], axis=1)
        del only_features_results_sorted_key2int
        del only_features_results_sorted

        #kdrew: remove self
        features_results_sorted_key2int_noself = features_results_sorted_key2int[features_results_sorted_key2int['key1_int'] != features_results_sorted_key2int['key2_int']]
        del features_results_sorted_key2int

        #kdrew: trim to just keys and prob
        keys_wprob = features_results_sorted_key2int_noself[['key1_int','key2_int','svm_pos_prob']]
        del features_results_sorted_key2int_noself

        #kdrew: remove duplicate rows where the keys are the same
        keys_nodups = keys_wprob.drop_duplicates(['key1_int','key2_int'])
        print keys_nodups
        print " "
        #kdrew: remove NAs
        keys_nodups = keys_nodups.dropna()
        print keys_nodups
        print " "
        #kdrew: convert type to int
        if args.int_convert:
            keys_nodups['key1_int'] = keys_nodups['key1_int'].astype(int)
            keys_nodups['key2_int'] = keys_nodups['key2_int'].astype(int)
        print keys_nodups
        print " "

    else:
        #kdrew: only want keys
        #keys_only = only_features_results_sorted_key2int_noself[['key1','key2']]
        keys_only = only_features_results_sorted_keys_noself[args.id_columns]
        #kdrew: remove duplicates
        keys_nodups = keys_only.drop_duplicates()
        print keys_nodups
        print " "
        #kdrew: remove NAs
        keys_nodups = keys_nodups.dropna()
        print keys_nodups
        print " "
        #kdrew: convert type to int
        if args.int_convert:
            keys_nodups[args.id_columns[0]] = keys_nodups[args.id_columns[0]].astype(int)
            #keys_nodups['key2_int'] = keys_nodups['key2_int'].astype(int)
            keys_nodups[args.id_columns[1]] = keys_nodups[args.id_columns[1]].astype(int)
        print keys_nodups
        print " "

    keys_nodups.to_csv(args.output_file,sep="\t",index=False,header=False)



if __name__ == "__main__":
    main()


