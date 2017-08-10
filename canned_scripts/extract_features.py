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
    elution = ef()
    elution.load(args.infile)
    feature_matrix = elution.extract_features(feature=args.feature,resampling=args.resampling,iterations=args.iterations,threshold=args.threshold)
    
    if feature_matrix is not None: # could be None if nothing left after thresholding, should do something else here though
        outfile = "_".join( [args.infile.split(".")[0], elution.analyses[ elution.analysis_count - 1 ]] ) + ".feat"
    
        feature_matrix.to_csv(outfile,index=False)
