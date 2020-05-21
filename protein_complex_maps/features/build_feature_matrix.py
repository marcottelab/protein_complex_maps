from __future__ import print_function
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import os.path

def build_matrix():

    parser = argparse.ArgumentParser(description="Builds feature matrix by concatenating pairs files (format: ID1 ID2 value")
    parser.add_argument("--input_pairs_files", action="store", dest="pairs_files", nargs='+', required=False, 
                                    help="Filenames of pairs files, space separated string")

    parser.add_argument("--input_pairs_list", action="store", dest="pairs_file_list", required=False, 
                                    help="File of pairs filenames, one per line")

    parser.add_argument("--sep", action="store", dest="sep", required=True, 
                                    help="Separator for input files")

    parser.add_argument("--prev_feature_matrix", action="store", nargs='?', dest='dfall', required=False)
    parser.add_argument("--store_interval", action="store", dest="interval", type=int, required=False, default=None, 
                                    help="If provided, stores an intermediate features matrix every x files joined")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--header", action="store", dest="header", required=False, default=True, 
                                    help="Do input pairs files have headers?")

    args = parser.parse_args()
    
    dfall = args.dfall

    if not args.pairs_files:
         if not args.pairs_file_list:
            print("either --input_pairs or --input_pairs_file required")
            return()

    print(dfall)
    if not dfall:
        dfall = pd.DataFrame(columns=['ID'])
     
    else:
        dfall = pd.read_csv(dfall)
        print(dfall)
        
    print("interval", args.interval) 

    #kdrew: read in pairs files
    dfall = dfall.set_index(['ID'])
    out_columns = []


    if args.pairs_files:
           pairs_files = args.pairs_files
    if args.pairs_file_list:
           pairs_files = pd.read_csv(args.pairs_file_list, header = None, names = ["filenames"])["filenames"].tolist()



    print(pairs_files)
    for i, filename in enumerate(pairs_files):
        print(i, filename)
        value_name = '.'.join(os.path.basename(filename).split('.')[:-2])
        out_columns.append(value_name)
     
        if args.header == True:
            df = pd.read_csv(filename, delimiter=args.sep, header=0, names=['A', 'B', value_name])
        elif args.header == False:
            df = pd.read_csv(filename, delimiter=args.sep, header=None, names=['A', 'B', value_name])


        print(df)
        df['ID']=df.apply(lambda x:'%s %s' % (x['A'],x['B']),axis=1)
        print("new column made")
        df_intermediate = df.drop(['A', 'B'], 1)
    
        df = df_intermediate
        print(df)

        if value_name in dfall.columns:
            dfall  = dfall.drop(value_name, 1)
        
        df = df.set_index(['ID'])
        print(df.shape)
        dfall_intermediate = dfall.join(df, how='outer', rsuffix=('_%s' % (i) ))
        dfall = None
        dfall = dfall_intermediate
        dfall_intermediate = None
        df = None
  
        # Save an intermediate file for long running joins
        if args.interval: 
            if i % args.interval == 0:
                if i > 0:
                    int_out_filename = args.out_filename +  ".intermediate_" + str(i) 
                    dfall.to_csv(int_out_filename)

           
        print(dfall.shape)

    dfall = dfall.reset_index()

    #kdrew: write to file 
    dfall.to_csv(args.out_filename, index=False)



if __name__ == "__main__":
    build_matrix()


