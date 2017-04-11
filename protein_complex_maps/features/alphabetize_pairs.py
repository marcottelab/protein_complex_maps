from __future__ import print_function
import argparse
import pandas as pd


#python ../../scripts/alphabetize_pairs.py --feature_pairs arathtraesorysjbraolselml_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob_c32_g0078125_pairs_noself_nodups_wprob.txt --outfile old_arathtraesorysjbraolselml_interactions.tmp

def alphabetized_check(df, column_ids, sample_size=1000):
    print(column_ids)
    df_sample = df.sample(sample_size)
    df_sample['alpha_id1'] = df_sample[column_ids].apply(sorted,axis=1)[column_ids[0]]
    ret_value = all(df_sample['alpha_id1'] == df_sample[column_ids[0]])
    return ret_value


def main():
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

    df.to_csv(args.outfile, header=False, index=False, sep=args.sep)





if __name__ == "__main__":
    main()


