from __future__ import print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

def build_matrix():

    parser = argparse.ArgumentParser(description="Builds feature matrix by concatenating pairs files")
    parser.add_argument("--input_pairs_files", action="store", dest="pairs_files", nargs='+', required=True, 
                                    help="Filenames of pairs files, format: id1\tid2\tvalue")
    parser.add_argument("--prev_feature_matrix", action="store", nargs='?', dest='dfall', required=False)
    parser.add_argument("--store_interval", action="store", dest="interval", type=int, required=False, default=None, 
                                    help="If provided, stores an intermediate features matrix every x files joined")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True, default=None, 
                                    help="Filename of output file")
    args = parser.parse_args()
    
    dfall = args.dfall



    print(dfall)
    if not dfall:
        dfall = pd.DataFrame(columns=['ID'])
     
    else:
   	 dfall = pd.read_table(dfall, sep=",")
         print(dfall)
        
         #dfall['pairset'] = map(frozenset, zip( dfall['ID1'].values, dfall['ID2'].values ))
         #a['IDHash'] = a.ID.apply(lambda r: tuple(sorted(r.iteritems())))
         #b['IDHash'] = b.ID.apply(lambda r: tuple(sorted(r.iteritems())))
         #print("dfall index set")
    #Claire: haven't checked this, but allow features to be added onto existing feature matrix. 
    #kdrew: read in pairs files
    #dfall = pd.DataFrame(columns=['pairset'])
    dfall = dfall.set_index(['ID'])
    out_columns = []


    for i, filename in enumerate(args.pairs_files):
        print(i, filename)
        value_name = '.'.join(os.path.basename(filename).split('.')[:-2])
        out_columns.append(value_name)
        df = pd.read_csv(filename, delimiter=' ', header=None, names=['A', 'B', value_name])
        print(df)
        df['ID']=df.apply(lambda x:'%s %s' % (x['A'],x['B']),axis=1)
        print("new column made")
        df_intermediate = df.drop(['A', 'B'], 1)
    
        df = df_intermediate
        print(df)

        #df.columns = ['ID',value_name])

        #print(filename)
        if value_name in dfall.columns:
            dfall  = dfall.drop(value_name, 1)
        
        #kdrew: create frozenset column with id1 and id2
        #df['pairset'] = map(frozenset, zip( df['id1'].values, df['id2'].values ))
        #print(df)
        #claire get rid of ID columns
        #df = df[['pairset', value_name]]
        df = df.set_index(['ID'])
        print(df.shape)

        dfall_intermediate = dfall.join(df, how='outer', rsuffix=('_%s' % (i) ))
        dfall = dfall_intermediate
        dfall_intermediate = None
        df = None
 
        if i % args.interval == 0:
             int_out_filename = "interval_" + str(i) + "_" + args.out_filename
             dfall.to_csv(int_out_filename)

           
        #print(dfall)
        #print(df)
        print(dfall.shape)
        #kdrew: outer merge dfs on frozenset column
        #dfall = dfall.merge(df, how='outer', on='pairset', suffixes=('', '_%s' % (i) ))

    dfall = dfall.reset_index()
    #dfall['ID1'] = [list(x)[0] for x in dfall['pairset'].values]
    #dfall['ID2'] = [list(x)[1] for x in dfall['pairset'].values]

    #kdrew: prepend global id columns to out vector
    #out_columns = ['ID1','ID2'] + out_columns

    #dfall = dfall.fillna(0.0)

    #kdrew: write to file 
    dfall.to_csv(args.out_filename, index=False)



if __name__ == "__main__":
    build_matrix()


