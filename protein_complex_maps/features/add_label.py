
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
    parser.add_argument("--input_negatives", action="store", dest="negatives", required=False, default=None,
                                    help="Filenanme of negative pairs, default = None (generated from processing positives)")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading csv, default=$")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="List of columns that specify ids in feature matrix (ex. gene_id bait_geneid")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--fillna", action="store", dest="fillna", required=False, default=None, 
                                    help="If set, fills NAs with input value")
    parser.add_argument("--int_convert", action="store_true", dest="int_convert", required=False, default=False, 
                                    help="Convert id_columns to int")
    parser.add_argument("--index_col0", action="store_true", dest="index_col0", required=False, default=False, 
                                    help="input_feature_matrix includes unnamed index column at position 0")
    args = parser.parse_args()


    if len(args.id_columns) != 2:
        print "Error: must provide two id columns"
        return -1

    if args.index_col0:
        feature_table = pd.read_csv(args.feature_matrix,sep=args.sep, index_col=0)
    else:
        feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)

    if args.fillna != None:
        feature_table = feature_table.fillna(float(args.fillna))

    positive_file = open(args.positives,"rb")

    all_proteins = set()
    ppis = set()
    neg_ppis = set()
    for line in positive_file.readlines():
        if len(line.split()) >= 2:
            id1 = line.split()[0]
            id2 = line.split()[1]
            all_proteins.add(id1)
            all_proteins.add(id2)

            ppi_str = str(sorted(list(frozenset([id1,id2]))))
            ppis.add(ppi_str)

    #kdrew: generate negative list by generating all pairs of proteins in positives but edges are not in positive list
    if args.negatives == None:
        for cpair in it.combinations(all_proteins,2):
            pair = str(sorted(list(frozenset(cpair))))
            if pair not in ppis:
                #print "pair is neg: %s" % ' '.join(pair)
                neg_ppis.add(pair)
    else:
        #kdrew: readin from file
        negative_file = open(args.negatives,"rb")
        for line in negative_file.readlines():
            if len(line.split()) >= 2:
                id1 = line.split()[0]
                id2 = line.split()[1]
                ppi_str = str(sorted(list(frozenset([id1,id2]))))
                neg_ppis.add(ppi_str)


    print "size of neg_ppis: %s" % len(neg_ppis)

    #kdrew: testing
    #ppis.add(frozenset([5987, 222389]))
    #neg_ppis.add(frozenset([10437, 64854]))

    #kdrew: the bioplex dataframe was weird when applying a set function because it would generate the set and put it in both columns as a dataframe
    #kdrew: it changed when the ids were floats instead of ints to a single column (not a dataframe)
    #is_ppis = feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id'].isin(ppis)
    #is_neg_ppis = feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id'].isin(neg_ppis)


    if 'id1_str' not in feature_table.columns and 'id2_str' not in feature_table.columns:
        if args.int_convert:
            feature_table['id1_str'] = feature_table[args.id_columns[0]].astype(int).apply(str)
            feature_table['id2_str'] = feature_table[args.id_columns[1]].astype(int).apply(str)
        else:
            feature_table['id1_str'] = feature_table[args.id_columns[0]].apply(str)
            feature_table['id2_str'] = feature_table[args.id_columns[1]].apply(str)
        feature_table['frozenset_ids'] = map(frozenset,feature_table[['id1_str','id2_str']].values)
        feature_table['frozenset_ids_str_order'] = feature_table['frozenset_ids'].apply(list).apply(sorted).apply(str)
    else:
        print "Warning: id1_str / id2_str are already in feature table"

    #print list(ppis)[:10]
    #print feature_table['frozenset_ids_str_order'].values[:10]

    #is_ppis = feature_table[args.id_columns].apply(set,axis=1).isin(ppis)
    #is_ppis = feature_table['frozenset_ids_order'].isin(ppis)
    is_ppis = [x in ppis for x in feature_table['frozenset_ids_str_order'].values] 
    #is_neg_ppis = feature_table[args.id_columns].apply(set,axis=1).isin(neg_ppis)
    #is_neg_ppis = feature_table['frozenset_ids_str_order'].isin(neg_ppis)
    is_neg_ppis = [x in neg_ppis for x in feature_table['frozenset_ids_str_order'].values]

    labels = [1 if is_ppis[index] else -1 if is_neg_ppis[index] else 0 for index in xrange(len(is_ppis))]
    
    #print feature_table[['gene_id','bait_geneid']].apply(set,axis=1)['gene_id']
    feature_table['label'] = labels
    feature_table['is_ppis'] = is_ppis
    feature_table['is_neg_ppis'] = is_neg_ppis
    print len(ppis)
    print len(neg_ppis)
    print is_ppis.count(True)
    print is_neg_ppis.count(True)
    #print feature_table[['gene_id','bait_geneid','is_ppis','is_neg_ppis','label']].head()


    #kdrew: weird extra column gets added, remove so it does not cause problems later
    if 'Unnamed: 0.1' in feature_table.columns:
        feature_table = feature_table.drop('Unnamed: 0.1',axis=1) 

    feature_table.to_csv(args.out_filename,sep=args.sep)

if __name__ == "__main__":
    main()


