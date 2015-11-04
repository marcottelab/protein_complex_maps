
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu


def main():

    parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--input_positives", action="store", dest="positives", required=True, 
                                    help="Filenanme of positive pairs")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading csv, default=$")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="List of columns that specify ids in feature matrix (ex. gene_id bait_geneid")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--fillna", action="store", dest="fillna", required=False, default=None, 
                                    help="If set, fills NAs with input value")
    args = parser.parse_args()


    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep, index_col=0)
    if args.fillna != None:
        feature_table = feature_table.fillna(float(args.fillna))

    positive_file = open(args.positives,"rb")

    all_proteins = set()
    ppis = set()
    neg_ppis = set()
    for line in positive_file.readlines():
        if len(line.split()) >= 2:
            id1 = int(line.split()[0])
            id2 = int(line.split()[1])
            all_proteins.add(id1)
            all_proteins.add(id2)

            ppis.add(frozenset([id1,id2]))

    for pair in it.combinations(all_proteins,2):
        if frozenset(pair) not in ppis:
            #print "pair is neg: %s" % ' '.join(pair)
            neg_ppis.add(frozenset(pair))


    print "size of neg_ppis: %s" % len(neg_ppis)

    #kdrew: testing
    #ppis.add(frozenset([5987, 222389]))
    #neg_ppis.add(frozenset([10437, 64854]))

    #kdrew: the bioplex dataframe was weird when applying a set function because it would generate the set and put it in both columns as a dataframe
    #kdrew: it changed when the ids were floats instead of ints to a single column (not a dataframe)
    #is_ppis = feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id'].isin(ppis)
    #is_neg_ppis = feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id'].isin(neg_ppis)
    is_ppis = feature_table[args.id_columns].apply(set,axis=1).isin(ppis)
    is_neg_ppis = feature_table[args.id_columns].apply(set,axis=1).isin(neg_ppis)

    labels = [1 if is_ppis[index] else -1 if is_neg_ppis[index] else 0 for index in xrange(len(is_ppis))]
    
    #print feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id']
    feature_table['label'] = labels
    feature_table['is_ppis'] = is_ppis
    feature_table['is_neg_ppis'] = is_neg_ppis
    print len(ppis)
    print len(neg_ppis)
    print is_ppis.sum()
    print is_neg_ppis.sum()
    #print feature_table[['gene_id','bait_geneid','is_ppis','is_neg_ppis','label']].head()


    #kdrew: weird extra column gets added, remove so it does not cause problems later
    if 'Unnamed: 0.1' in feature_table.columns:
        feature_table = feature_table.drop('Unnamed: 0.1',axis=1) 

    feature_table.to_csv(args.out_filename,sep=args.sep)

if __name__ == "__main__":
    main()


