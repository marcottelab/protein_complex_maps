from __future__ import print_function
import sys
import ast
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.features.alphabetize_pairs as ap

def load_ppis(ppi_file, ppi_sep, id_sep, label=1):

    ppis = pd.DataFrame(pd.read_table(ppi_file, sep=ppi_sep, header=None))
    ppis.columns = ['ID1','ID2']
    if not ap.alphabetized_check(ppis, ['ID1','ID2'], 100):
        sys.exit("ERROR: input_positives are not alphabetized, please run alphabetize_pairs.py")
    ppis['label'] = label
    ppis['ID'] = ppis['ID1'] + id_sep + ppis['ID2']
    return ppis

def main():

    parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--input_positives", action="store", dest="positives", required=True, 
                                    help="Filename of positive pairs")
    parser.add_argument("--input_negatives", action="store", dest="negatives", required=False, default=None,
                                    help="Filename of negative pairs, default = None (generated from processing positives)")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading csv, default=$")
    parser.add_argument("--ppi_sep", action="store", dest="ppi_sep", required=False, default=',',
                                    help="Separator for input_positives and input_negatives files, default=,")
    parser.add_argument("--id_column", action="store", nargs='+', dest="id_column", required=False, default = ["ID"], 
                                    help="column(s) that specifies the id in feature matrix (ex. ID), can pass in multiple columns")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--fillna", action="store", dest="fillna", required=False, default=None, 
                                    help="If set, fills NAs with input value")
    parser.add_argument("--id_sep", action="store", dest="id_sep", type=str, required=False, default="_", 
                                    help="If ID is stored in one column, give the spacing used between pairs ie. id1_id2")
    parser.add_argument("--sel_columns", action="store", dest="cols", type=str, required=False, default=None, 
                                    help="Optional columns to save")   
    #parser.add_argument("--int_convert", action="store_true", dest="int_convert", required=False, default=False, 
    #                                help="Convert id_column to int")
    #parser.add_argument("--index_col0", action="store_true", dest="index_col0", required=False, default=False, 
    #                                help="input_feature_matrix includes unnamed index column at position 0")
    args = parser.parse_args()

    pos_ppis = load_ppis(args.positives, args.ppi_sep, args.id_sep, 1)
    print("size of pos_ppis: %s" % len(pos_ppis))
    if args.negatives:
        neg_ppis = load_ppis(args.negatives, args.ppi_sep, args.id_sep, -1)
        print("size of neg_ppis: %s" % len(neg_ppis))
        all_ppis = pd.DataFrame(pd.concat([pos_ppis, neg_ppis]))

    else:
        all_ppis = pos_ppis

    print(all_ppis)
    print(type(all_ppis))

    all_ppis = all_ppis.set_index(['ID'])
    print(all_ppis.head)
  
    feature_table = pd.DataFrame(pd.read_table(args.feature_matrix, sep=args.sep))

    if not ap.alphabetized_check(feature_table, args.id_column):
        sys.exit("ERROR: feature_matrix is not alphabetized, please run alphabetize_pairs.py")
    #kdrew: if multiple ids were passed in, combine into a string
    feature_table['ID'] = feature_table[args.id_column].apply(lambda x: args.id_sep.join(map(str,x)), axis=1)

    feature_table = feature_table.set_index(['ID'])

    print(feature_table.head)
    print(type(feature_table))    

    if args.fillna != None:
        feature_table = feature_table.fillna(float(args.fillna))



    
    labeled_feature_table = feature_table.join(all_ppis, how="left", rsuffix="_x")   

    print("pos/neg")

    labeled_feature_table['label'] = labeled_feature_table['label'].fillna(0)


    #kdrew: weird extra column gets added, remove so it does not cause problems later
    #CDM probably unnecessary now
    if 'Unnamed: 0.1' in labeled_feature_table.columns:
        labeled_feature_table = labeled_feature_table.drop('Unnamed: 0.1',axis=1) 


    labeled_feature_table = labeled_feature_table.reset_index()
    if args.cols:
        print(args.cols)
        collist = args.cols.split(' ')
        labeled_feature_table = labeled_feature_table[collist]


    labeled_feature_table.to_csv(args.out_filename,sep=args.sep, index=False)

if __name__ == "__main__":
    main()


