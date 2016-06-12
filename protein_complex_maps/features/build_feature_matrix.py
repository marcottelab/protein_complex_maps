
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

def main():

    parser = argparse.ArgumentParser(description="Builds feature matrix by concatenating pairs files")
    parser.add_argument("--input_pairs_files", action="store", dest="pairs_files", nargs='+', required=True, 
                                    help="Filenames of pairs files, format: id1\tid2\tvalue")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True, default=None, 
                                    help="Filename of output file")
    args = parser.parse_args()



    #kdrew: read in pairs files
    dfall = pd.DataFrame(columns=['pairset'])
    out_columns = []


    for i, filename in enumerate(args.pairs_files):

        value_name = '.'.join(os.path.basename(filename).split('.')[:-1])
        out_columns.append(value_name)
        df = pd.read_csv(filename, delimiter='\t', header=None, names=['id1','id2',value_name])
        
        #kdrew: create frozenset column with id1 and id2
        df['pairset'] = map(frozenset, zip( df['id1'].values, df['id2'].values ))

        #kdrew: outer merge dfs on frozenset column
        dfall = dfall.merge(df, how='outer', on='pairset', suffixes=('', '_%s' % (i) ))


    dfall['ID1'] = [list(x)[0] for x in dfall['pairset'].values]
    dfall['ID2'] = [list(x)[1] for x in dfall['pairset'].values]

    #kdrew: prepend global id columns to out vector
    out_columns = ['ID1','ID2'] + out_columns

    #dfall = dfall.fillna(0.0)

    #kdrew: write to file 
    dfall.to_csv( args.out_filename, columns=out_columns, index=False)



if __name__ == "__main__":
    main()


