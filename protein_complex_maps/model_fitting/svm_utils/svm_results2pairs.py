
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

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
    parser.add_argument("--nofilter_label", action="store_true", dest="nofilter_label", required=False, default=False,
                                    help="Do not filter by label, combine the entire feature matrix with results, default=False")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=False, default=['key1','key2'],
                                    help="List of columns that specify ids in feature matrix, default: key1 key2")
    parser.add_argument("--id_sep", action="store", dest="id_sep", required=False, default=' ',
                                    help="If IDs are stored in one column, separator to split on")

    parser.add_argument("--int_convert", action="store_true", dest="int_convert", required=False, default=False, 
                                    help="Convert id_columns to int, default=False")
    parser.add_argument("--add_label", action="store_true", dest="add_label", required=False, default=False, 
                                    help="Add the original training label, default=False")


    parser.add_argument("--label_column", action="store", dest="label_column",  required=False, default='label',
                                    help="Column containing label")


    args = parser.parse_args()

    results = pd.read_csv(args.input_results,sep=' ')
    results.columns = ['svm_predicted_label','svm_pos_prob','svm_neg_prob']    


    print results
    print " "

    features = pd.read_csv(args.feature_matrix,sep=args.sep)
    print features
    print " "
    if args.nofilter_label:
        only_features = features
    else:
        if args.label_not0:
            only_features = features[features['label'] != 0]
        else:
            only_features = features[features['label'] == 0]
    only_features.index = range(len(only_features))

    #cmcwhite keep only id and label from feature matrix
    if args.label_column:
        only_feature_ids = only_features[args.id_columns + [args.label_column]]  
    else:
        only_feature_ids = only_features[args.id_columns]  

    key1 = args.id_columns[0]
    key2 = args.id_columns[1]
      

    print only_features
    print " "

    print len(results)
    print len(only_feature_ids)
    assert len(results) == len(only_feature_ids), "features and results are not same length"
    #kdrew: put the results and the features together
    print("Begin concat")
    only_features_results = pd.concat([results,only_feature_ids],axis=1)
    print only_features_results
    print " "

    #kdrew: cleaning up space
    del results
    del only_features

    only_features_results = only_features_results.astype('object')
    print only_features_results
    print " "

    print("sorting")
    only_features_results_sorted = only_features_results.sort_values(by='svm_pos_prob',ascending=False)
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
    only_features_results_sorted_keys_noself = only_features_results_sorted_keys[only_features_results_sorted_keys[key1] != only_features_results_sorted_keys[key2]]

    legacy = False
    if args.add_prob:
        if legacy:
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

            if args.add_label:
                keys_and_prob = only_features_results_sorted[args.id_columns+[args.label_column]+['svm_pos_prob']]

            else:
                keys_and_prob = only_features_results_sorted[args.id_columns+['svm_pos_prob']]

            keys_and_prob_noself = keys_and_prob[keys_and_prob[key1] != keys_and_prob[key2]]
            keys_nodups = keys_and_prob_noself.drop_duplicates()
            keys_nodups = keys_nodups.dropna()


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
            keys_nodups[key1] = keys_nodups[key1].astype(int)
            #keys_nodups['key2_int'] = keys_nodups['key2_int'].astype(int)
            keys_nodups[key2] = keys_nodups[key2].astype(int)
        print keys_nodups
        print " "

    keys_nodups.to_csv(args.output_file,sep="\t",index=False,header=False)



if __name__ == "__main__":
    main()


