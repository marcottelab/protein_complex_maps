#! /usr/bin/env python

import argparse
from protein_complex_maps.features.ExtractFeatures.Features import ElutFeatures as ef

parser = argparse.ArgumentParser(description="Extract features from a fractionation mass-spec experiment")
parser.add_argument("infile")
parser.add_argument("--format", default='tsv', choices=['tsv','csv'], help="Format of infile: 'tsv' or 'csv'")
parser.add_argument("-f", "--feature", default="pearsonR",choices=ef.available_features, help="Feature to extract")
parser.add_argument("-r", "--resampling", default=False, choices=ef.resampling_strategies, help="Resample the data and average the feature")
parser.add_argument("-i", "--iterations", default=False, type=int, help="Number of iterations for the resampling strategy")
parser.add_argument("-t", "--threshold", default=False, type=int, help="Remove proteins with fewer than this number of PSMs")
parser.add_argument("-n", "--normalize", nargs="+", default=False, choices=["row_sum", "row_max","column"], help="Whether to normalize by row or column")
parser.add_argument("--as_pickle", action='store_true', help='Save output as pickled DataFrame')
args = parser.parse_args()

if __name__ == '__main__':
    elution = ef()
    elution.load(args.infile, format=args.format)
    feature_matrix = elution.extract_features(feature=args.feature,resampling=args.resampling,iterations=args.iterations,threshold=args.threshold,normalize=args.normalize)
    
    if feature_matrix is not None: # could be None if nothing left after thresholding, should do something else here though
        outfile = "_".join( [args.infile.split(".")[0], elution.analyses[ elution.analysis_count - 1 ]] ) + ".feat"
    
        if args.as_pickle:
            feature_matrix.to_pickle(outfile + ".p")
        else:
            feature_matrix.to_csv(outfile,index=False)
