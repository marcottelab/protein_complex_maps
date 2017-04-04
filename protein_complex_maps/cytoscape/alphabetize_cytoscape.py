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


    parser = argparse.ArgumentParser(description="Alphabetize")
    parser.add_argument("--interaction_network", action="store", dest= "df", required=True)
    parser.add_argument("--outfile", action="store", required=True)
    args = parser.parse_args()
    df = pd.read_table(args.df, sep='\t')

    df.columns=['ID1_num', 'ID2_num', 'ID1', 'ID2']


    df['ID'] = map(sorted, zip(df['ID1'].values, df['ID2'].values))

    df = df.drop(['ID1', 'ID2'], axis=1)

   
    #Change list format entries to string entries 
    df['ID'] = df['ID'].apply(lambda x: ' '.join(x))
  


    df.to_csv(args.outfile, index=False)



if __name__ == "__main__":
    order_identifiers()


