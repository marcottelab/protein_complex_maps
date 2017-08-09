#! /usr/bin/env python

import argparse
from elut import ElutFeatures as ef

parser = argparse.ArgumentParser(description="Extract features from a fractionation mass-spec experiment")
parser.add_argument("infile")
parser.add_argument("-f", "--feature", default="pearsonR",choices=ef.available_features)
parser.add_argument("-r", "--resampling", default=None, choices=ef.resampling_strategies)
parser.add_argument("-i", "--iterations", default=None, type=int)
parser.add_argument("-t", "--threshold", default=None, type=int)
args = parser.parse_args()

if __name__ == '__main__':
    ## This won't work if threshold or iterations are None. Fix
    outfile = "_".join([args.infile.split(".")[0], args.feature, args.resampling, str(args.iterations)+"reps", "thresh"+str(args.threshold)]) + ".feat"
    elution = ef()
    elution.load(args.infile)
    feature_matrix = elution.extract_features(feature=args.feature,resampling=args.resampling,iterations=args.iterations,threshold=args.threshold)
    feature_matrix.to_csv(outfile,index=False)
