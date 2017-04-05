from __future__ import print_function
import argparse
import pandas as pd


#python ../../scripts/alphabetize_pairs.py --feature_pairs arathtraesorysjbraolselml_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob_c32_g0078125_pairs_noself_nodups_wprob.txt --outfile old_arathtraesorysjbraolselml_interactions.tmp


def order_identifiers():
    '''
    This function takes a file of pairwise ID features:
    A B value
    D C value
    And orders the identifiers:
    A B value
    C D value
    This is because different feature makers may order the pairwise IDs differently
    '''


    parser = argparse.ArgumentParser(description="Selects multiple selections of columns from a list")
    parser.add_argument("--feature_pairs", action="store", dest= "df", required=True)
    parser.add_argument("--outfile", action="store", required=True, 
                                    help="Filename of comma separated output")
    parser.add_argument("--sep", action = "store",dest='sep' , required=True)
    parser.add_argument("--columns", nargs='+', action="store", type=int, dest="columns2alphabetize", required=False, default=[0,1], 
                                    help="List of columns to alphabetize, default = 0,1")
    args = parser.parse_args()
    df = pd.read_table(args.df, sep=args.sep, header=None)
    df[df.columns[args.columns2alphabetize]] = df[df.columns[args.columns2alphabetize]].apply(sorted,axis=1)

    df.to_csv(args.outfile, header=False, index=False)



if __name__ == "__main__":
    order_identifiers()


