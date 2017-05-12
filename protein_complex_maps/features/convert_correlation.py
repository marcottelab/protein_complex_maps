
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

def main():

    parser = argparse.ArgumentParser(description="Converts correlation matrix from elution profile to pairwise list")
    parser.add_argument("--input_correlation_matrix", action="store", dest="correlation_matrix", required=True, 
                                    help="Filename of correlation matrix, ex. output of blake_complexes/score.py")
    parser.add_argument("--input_elution_profile", action="store", dest="elution_profile", required=True, 
                                    help="Filename of elution profile for protein ids")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=True, default=None, 
                                    help="Filename of output file")
    parser.add_argument("--noheader", action="store_true", dest="noheader", required=False, default=False, 
                                    help="Flag to set if elution profile does not have a header, default=False")
    parser.add_argument("--sep", action="store_true", dest="sep", required=False, default=' ', 
                                    help="Separator for the elution profile")


    args = parser.parse_args()


    protein_ids = []
    correlation_scores = dict()
    sep=args.sep
    ep_file = open(args.elution_profile,"rb")
    for i, line in enumerate(ep_file.readlines()):
        if i > 0 or args.noheader:
            protein_ids.append(line.split(sep)[0].strip())

    ep_file.close()


    #print protein_ids

    cm_file = open(args.correlation_matrix,"rb")
    for i, line in enumerate(cm_file.readlines()):
        for j, val in enumerate(line.split()):
            correlation_scores[frozenset([protein_ids[i],protein_ids[j]])] = float(val)

    cm_file.close()


    sorted_scores = sorted(correlation_scores.items(), key=operator.itemgetter(1), reverse=True)

    out_file = open(args.out_filename, "wb")
    for s, score in sorted_scores:

        #kdrew: don't output pairs where the ids are the same
        if len(list(s)) != 2:
            continue

        out_file.write("%s\t%s\t%s\n" % (list(s)[0], list(s)[1], score))

    out_file.close()



if __name__ == "__main__":
    main()


