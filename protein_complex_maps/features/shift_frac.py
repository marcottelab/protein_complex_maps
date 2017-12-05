
import argparse
import numpy as np
import pandas as pd
import protein_complex_maps.features.ExtractFeatures.Features as eff


def main():

    parser = argparse.ArgumentParser(description="Calculate shift score between two fractionation experiments")
    parser.add_argument("--elution_files", action="store", nargs='+', dest="elution_files", required=True, 
                                    help="Elution files (.elut)")
    parser.add_argument("--normalize", action="store_true", dest="normalize", required=False, default=False,
                                    help="Normalize sum by the total number of counts in both elutions")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")

    args = parser.parse_args()
    
    elutions = []
    for efile in args.elution_files:
        elut = eff.Elut()
        elut.load(efile,format='tsv')
        elutions.append(elut)

    if len(elutions) >= 2:
        shift_frac_scores = calc_shift_frac(elutions[0], elutions[1], normalize_totalCounts=args.normalize)
        print shift_frac_scores
        if args.out_filename != None:
            shift_frac_scores.sort_values(ascending=False).to_csv(args.out_filename)

def calc_shift_frac(elut1, elut2, normalize_totalCounts=False):

    #kdrew: set columns to be the same, might want to do some error checking to ensure lengths match, also if any realignment is necessary this is the place to do it.
    elut2.df.columns = elut1.df.columns

    #kdrew: add empty rows for the ids in elut1 that are not in elut2 and vice versa
    elut1_ids = set(elut1.df.index)
    elut2_ids = set(elut2.df.index)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut1_not_elut2_ids = elut1_ids - elut2_ids
    elut1_not_elut2 = elut1.df.loc[list(elut1_not_elut2_ids)]
    elut1_not_elut2[:] = 0.0
    elut2.df = elut2.df.append(elut1_not_elut2)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut2_not_elut1_ids = elut2_ids - elut1_ids
    elut2_not_elut1 = elut2.df.loc[list(elut2_not_elut1_ids)]
    elut2_not_elut1[:] = 0.0
    elut1.df = elut1.df.append(elut2_not_elut1)

    elut_diff = elut1.df.subtract(elut2.df)
    shift_frac_sum = np.abs(elut_diff.sum(axis='columns'))

    if normalize_totalCounts:
        shift_frac_sum = shift_frac_sum * shift_frac_sum/(elut1.df.sum(axis='columns') + elut2.df.sum(axis='columns'))

    return shift_frac_sum

if __name__ == "__main__":
    main()

