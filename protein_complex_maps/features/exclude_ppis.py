
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it

import protein_complex_maps.correlation_util as cu
import protein_complex_maps.protein_util as pu


def main():

    parser = argparse.ArgumentParser(description="Removes specified ppis from feature matrix")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading csv, default=$")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=True, 
                                    help="List of columns that specify ids in feature matrix (ex. gene_id bait_geneid")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--exclude_file", action="store", dest="exclude_file", required=True, 
                                    help="File of protein pairs which should be excluded from output")

    args = parser.parse_args()


    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)

    exclude_file = open(args.exclude_file,"rb")

    ppis = set()
    for line in exclude_file.readlines():
        if len(line.split()) >= 2:
            id1 = int(line.split()[0])
            id2 = int(line.split()[1])
            ppis.add(frozenset([id1,id2]))

    is_ppis = feature_table[args.id_columns].apply(set,axis=1).isin(ppis)
    print is_ppis

    feature_table = feature_table.drop(feature_table[is_ppis].index)

    #kdrew: weird extra column gets added, remove so it does not cause problems later
    if 'Unnamed: 0.1' in feature_table.columns:
        feature_table = feature_table.drop('Unnamed: 0.1',axis=1) 

    feature_table.to_csv(args.out_filename,sep=args.sep)

if __name__ == "__main__":
    main()


