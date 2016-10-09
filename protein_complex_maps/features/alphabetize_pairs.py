from __future__ import print_function
import argparse
import pandas as pd


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


    parser = argparse.ArgumentParser(description="Selects multiiple selections of columns from a list")
    parser.add_argument("--feature_pairs", action="store", dest= "df", required=True)
    parser.add_argument("--outfile", action="store", required=True)
    args = parser.parse_args()
 
    print("opening df")
    df = pd.read_table(args.df, sep="\t", header=None)
    print(df)
    df.columns = ['A', 'B', 'corr']

    print(df)

    df['ID'] = map(sorted, zip(df['A'].values, df['B'].values))

    df = df.drop(['A', 'B'], axis=1)


    df = df[['ID', 'corr']]

    

    df.to_csv(args.outfile, sep='\t', header=False, index=False)


    print(df)

if __name__ == "__main__":
    order_identifiers()


