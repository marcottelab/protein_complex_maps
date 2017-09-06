from __future__ import print_function
import argparse
import pandas as pd
#sys.path.append('/project/cmcwhite/protein_complex_maps/protein_complex_maps')
#import protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )

    parser.add_argument("--pairs_filename", action="store", type=str, dest="pairs_filename", required=True,
                                            help="PairWclustID")


    parser.add_argument("--train_filename", action="store", dest="train_filename", required=True,
                                            help="Filename of corum train interactions with label, csv.")
    parser.add_argument("--test_filename", action="store", dest="test_filename", required=True,
                                            help="Filename of corum test interactions with label, csv.")


    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")
    args = parser.parse_args()

    pairs_df = pd.read_csv(args.pairs_filename, index_col=False, sep=",")
  
    train_df = pd.read_csv(args.train_filename, index_col=False, sep="\t", header=None)
    #print(train_df)

    train_df.columns = ['ID', 'train_class']

    test_df = pd.read_csv(args.test_filename, index_col=False, sep="\t", header=None)
    test_df.columns = ['ID','test_class']

    print(pairs_df)


    pairs_df = pairs_df.set_index(['ID'])
    train_df = train_df.set_index(['ID'])
    test_df = test_df.set_index(['ID'])

    annotpairs_df = pairs_df.join(train_df, how="left")
    annotpairs2_df = annotpairs_df.join(test_df, how="left")


    annotpairs2_df = annotpairs2_df.reset_index()

    annotpairs2_df.to_csv(args.output_filename, index=False, sep="\t")
if __name__ == "__main__":
    main()
