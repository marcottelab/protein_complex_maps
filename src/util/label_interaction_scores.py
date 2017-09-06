from __future__ import print_function
import argparse
import pandas as pd
#sys.path.append('/project/cmcwhite/protein_complex_maps/protein_complex_maps')
#import protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )

    parser.add_argument("--pairs_filename", action="store", type=str, dest="pairs_filename", required=True,
                                            help="orderedpairsWclustIDtraintest")


    parser.add_argument("--scores_filename", action="store", dest="scores_filename", required=True,
                                            help="Filename of pairs and svm scores . Filtered05.txt")


    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")
    args = parser.parse_args()

    pairs_df = pd.read_csv(args.pairs_filename, index_col=False, sep="\t")
  
    scores_df = pd.read_csv(args.scores_filename, index_col=False, sep=" ", header=None)
    print(scores_df)

    scores_df.columns = ['ID1', 'ID2', 'svmscores']

    scores_df['ID'] = map(sorted, zip(scores_df['ID1'].values, scores_df['ID2'].values))

    scores_df = scores_df.drop(['ID1', 'ID2'], axis=1)

    scores_df['ID'] = scores_df['ID'].apply(lambda x: ' '.join(x))


    pairs_df = pairs_df.set_index(["ID"])
 
    

    scores_df = scores_df.set_index(["ID"])




    scorepairs_df = pairs_df.join(scores_df, how="left")


    scorepairs_df = scorepairs_df.reset_index()

    scorepairs_df.to_csv(args.output_filename, index=False, sep="\t")
if __name__ == "__main__":
    main()
