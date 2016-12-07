from __future__ import print_function
import argparse
import numpy as np
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
    parser.add_argument("--input_pairs", action="store", dest="pairs", required=True,
                                    help="Filenanme of positive pairs")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None,
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--label", action="store", dest="label", required=False, default=None,
                                    help="Fills a new column with a value")
    parser.add_argument("--sep", action="store", dest="sep", required=True, default=',',
                                    help="sep for the infile and outfile")
 
    args = parser.parse_args()


    #if len(args.id_column) != 2:
    #    print "Error: must provide two id columns"
    #    return -1

    #if args.index_col0:
    #    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep, index_col=0)
    #else:
    #    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)


    ppis = pd.DataFrame(pd.read_table(args.pairs, sep=args.sep, header=None))
    ppis.columns = ['A', 'B']
    print(ppis) 
    ppis['ID'] = map(sorted, zip(ppis['A'].values, ppis['B'].values))
    ppis = ppis.drop(['A', 'B'], axis=1) 
    ppis = ppis[['ID']]
    print(ppis)
    ppis['ID']= ppis['ID'].astype(str)
    print(ppis)
    ppis['ID'] = ppis['ID'].str.replace('[', '')
    ppis['ID'] = ppis['ID'].str.replace(']', '')
    ppis['ID'] = ppis['ID'].str.replace(',', '')
    ppis['ID'] = ppis['ID'].str.replace("'", "")
   

    print("size of pos_ppis: %s" % len(ppis))

    if args.label:

        ppis['label'] = args.label

    print(ppis)

    if args.out_filename:

        ppis.to_csv(args.out_filename, sep='\t', index=False)


main()
