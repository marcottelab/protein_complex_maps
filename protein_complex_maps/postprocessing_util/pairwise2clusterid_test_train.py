

import argparse
import itertools as it
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Tool to create table of test and training PPIs to pairwise interactions (and scores) with corresponding cluster_ids" )
    parser.add_argument("--pairwise_filename", action="store", dest="pairwise_filename", required=True,
                                            help="Filename of pairwise file, (ie. one line per ppi, third column contains score)")
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--train_ppis_filename", action="store", dest="train_ppis_filename", required=True,
                                            help="Filename of train ppis file, (ie. one line per ppi)")
    parser.add_argument("--test_ppis_filename", action="store", dest="test_ppis_filename", required=True,
                                            help="Filename of test ppis file, (ie. one line per ppi)")
    parser.add_argument("--train_neg_ppis_filename", action="store", dest="train_neg_ppis_filename", required=True,
                                            help="Filename of train negative ppis file, (ie. one line per ppi)")
    parser.add_argument("--test_neg_ppis_filename", action="store", dest="test_neg_ppis_filename", required=True,
                                            help="Filename of test negative ppis file, (ie. one line per ppi)")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True,
                                    help="Output filename ")
    args = parser.parse_args()

    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())

    pairwise_interactions = dict()
    f2 = open(args.pairwise_filename, "rb")
    for line in f2.readlines():
        try:
            pairwise_interactions[frozenset([line.split()[0],line.split()[1]])] = line.split()[2]
        except IndexError:
            continue

    test_positive_interactions = set()
    f3 = open(args.test_ppis_filename, "rb")
    for line in f3.readlines():
        try:
            test_positive_interactions.add(frozenset([line.split()[0],line.split()[1]]))
        except IndexError:
            continue

    train_positive_interactions = set()
    f4 = open(args.train_ppis_filename, "rb")
    for line in f4.readlines():
        try:
            train_positive_interactions.add(frozenset([line.split()[0],line.split()[1]]))
        except IndexError:
            continue

    test_negative_interactions = set()
    f5 = open(args.test_neg_ppis_filename, "rb")
    for line in f5.readlines():
        try:
            test_negative_interactions.add(frozenset([line.split()[0],line.split()[1]]))
        except IndexError:
            continue

    train_negative_interactions = set()
    f6 = open(args.train_neg_ppis_filename, "rb")
    for line in f6.readlines():
        try:
            train_negative_interactions.add(frozenset([line.split()[0],line.split()[1]]))
        except IndexError:
            continue

    print "read input files"
    fout = open(args.out_filename,"wb")

    fout.write("id1\tgoldstandard\n")

    print "Wrote Header"

    for clustid, cluster in enumerate(clusters):
        for prot_pair in it.combinations(cluster,2):
            goldstandard = ""
            id1 = "%s_%s" % (clustid, prot_pair[0])
            id2 = "%s_%s" % (clustid, prot_pair[1])
            try:
                if frozenset([prot_pair[0],prot_pair[1]]) in test_positive_interactions:
                    goldstandard = "test_pos"
                elif frozenset([prot_pair[0],prot_pair[1]]) in train_positive_interactions:
                    goldstandard = "train_pos"
                elif frozenset([prot_pair[0],prot_pair[1]]) in train_negative_interactions:
                    goldstandard = "train_neg"
                elif frozenset([prot_pair[0],prot_pair[1]]) in test_negative_interactions:
                    goldstandard = "test_neg"
            except KeyError:
                continue

            fout.write("%s (pp) %s\t%s" % (id1, id2, goldstandard))
            fout.write("\n")

    fout.close()

if __name__ == "__main__":
    main()
