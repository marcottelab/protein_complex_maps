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
    parser.add_argument("--outfile", action="store", required=True)
    parser.add_argument("--sep", action = "store",dest='sep' , required=True)
    parser.add_argument("--numcol", action="store", type=int, dest="valcol", required=False, default=2, help="If there is a third column beyond pairwise ids")
    args = parser.parse_args()
    print(args.valcol)
    if args.valcol == 3:
        print("Three columns")


    print("opening df")
    df = pd.read_table(args.df, sep=args.sep, header=None)
    print(df)
    if args.valcol == 3:

        df.columns = ['A', 'B', 'corr']

    else:
        df.columns = ['A', 'B']

    print(df)

    df['ID'] = map(sorted, zip(df['A'].values, df['B'].values))

    df = df.drop(['A', 'B'], axis=1)

    if args.valcol==3:
        df = df[['ID', 'corr']]

    else:
        df = df[['ID']]
   
    #Change list format entries to string entries 
    df['ID'] = df['ID'].apply(lambda x: ' '.join(x))
  


    df.to_csv(args.outfile, header=False, index=False)



if __name__ == "__main__":
    order_identifiers()


