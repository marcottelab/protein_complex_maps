from __future__ import print_function
import sys
import ast
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.features.alphabetize_pairs as ap


def load_ppis(positives, negatives, ppi_sep, id_sep):

    pos_ppis = pd.DataFrame(pd.read_table(positives, sep=ppi_sep, header=None, engine="python", converters={0: str, 1: str}))
    pos_ppis.columns = ['ID1','ID2']
    print(pos_ppis)

    #if args.check_alphabetized:
    #   if not ap.alphabetized_check(pos_ppis, ['ID1','ID2']):
    #     print("alphabetizing positives ppis")
    #     pos_ppis = ap.alphabetize_df(pos_ppis, [list(pos_ppis.columns).index(x) for x in ['ID1','ID2']])
    #     if not ap.alphabetized_check(pos_ppis, ['ID1','ID2']):
    #         sys.exit("ERROR: input_positives are not alphabetized, please run alphabetize_pairs.py")

    print(pos_ppis)
    print(pos_ppis['ID1'])
    print(pos_ppis['ID2'])
    print("size of pos_pos_ppis: %s" % len(pos_ppis))
    pos_ppis['label'] = 1
    pos_ppis['ID'] = pos_ppis['ID1'] + id_sep + pos_ppis['ID2']

    neg_ppis = pd.DataFrame(pd.read_table(negatives, sep=ppi_sep, header=None, engine="python", converters={0: str, 1: str}))
    neg_ppis.columns = ['ID1','ID2']
    #if not ap.alphabetized_check(neg_ppis, ['ID1','ID2'], 10):
    #print("alphabetizing negatives ppis")
    #    neg_ppis = ap.alphabetize_df(neg_ppis, [list(neg_ppis.columns).index(x) for x in ['ID1','ID2']])
       # if not ap.alphabetized_check(neg_ppis, ['ID1','ID2']):
       #     sys.exit("ERROR: input_negatives are not alphabetized, please run alphabetize_pairs.py")
    print(neg_ppis)
    neg_ppis['label'] = -1
    neg_ppis['ID'] = neg_ppis['ID1'] + id_sep + neg_ppis['ID2']
   
    
    print("size of neg_ppis: %s" % len(neg_ppis))

    
    all_ppis = pd.DataFrame(pd.concat([pos_ppis, neg_ppis]))
    all_ppis = all_ppis.drop(['ID1', 'ID2'], axis=1)
    print(all_ppis)
    print(type(all_ppis))

    all_ppis = all_ppis.set_index(['ID'])
    print(all_ppis.head)
    return(all_ppis)

def label_table(feature_table, id_column, id_sep, all_ppis, fillna, cols):

    #kdrew: if multiple ids were passed in, combine into a string
    feature_table['ID'] = feature_table[id_column].apply(lambda x: id_sep.join(map(str,x)), axis=1)

    feature_table = feature_table.set_index(['ID'])

    if cols:
        print(cols)
        collist = cols.split(' ')
        labeled_feature_table = labeled_feature_table[collist]

    if fillna != None:
        feature_table = feature_table.fillna(float(fillna))


    labeled_feature_table = feature_table.join(all_ppis, how="left", rsuffix="_x")   
    labeled_feature_table['label'] = labeled_feature_table['label'].fillna(0)

    return(labeled_feature_table)

 

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
    parser.add_argument("--check_alphabetized", action="store_true", required=False, default=False, 
                                    help="Whether to perform check that input IDs are alphabetized")
    parser.add_argument("--chunksize", action="store", dest="chunksize", type=int, required=False, default=100000,
                                    help="Chunksize for processing, default = 100000")    
   
    args = parser.parse_args()

    all_ppis = load_ppis(args.positives, args.negatives, args.ppi_sep, args.id_sep)

    #iterator = pd.read_csv(args.df, sep=args.sep, header=0, engine="python",  chunksize=args.chunksize, iterator = True)
    #feature_table = pd.DataFrame(pd.read_table(args.feature_matrix, sep=args.sep, engine="python"))
    iterator = pd.read_table(args.feature_matrix, sep=args.sep, engine="python", chunksize = args.chunksize, iterator = True)


    firstpass = True
    for feature_table in iterator:
        labeled_feature_table(feature_table, args.id_column, args.id_sep, all_ppis, args.fillna, args.cols)

        if firstpass == True:
            #alphabetized_df.to_csv(args.outfile, header=True, index=False, sep=args.sep)
            labeled_feature_table.to_csv(args.out_filename,sep=args.sep, index=True, header = True)
        else: # append without writing the header
            labeled_feature_table.to_csv(args.out_filename,sep=args.sep, index=True, header = False, mode = 'a')
        
        firstpass = False
 
if __name__ == "__main__":
    main()


