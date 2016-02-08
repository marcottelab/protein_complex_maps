

import argparse
import itertools as it
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Tool to create table to pairwise interactions (and scores) with corresponding cluster_ids" )
    parser.add_argument("--pairwise_filename", action="store", dest="pairwise_filename", required=True,
                                            help="Filename of pairwise file, (ie. one line per ppi, third column contains score)")
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=False, default=None, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=',',
                                    help="Column separator for input feature matrix file, default=,")
    parser.add_argument("--header_names", action="store", nargs='+', dest="header_names", required=False, default=None,
                                    help="Names for header relating to features")
    parser.add_argument("--features", action="store", nargs='+', dest="features", required=False, default=None,
                                    help="Names of fields for features in input_feature_matrix, will test !=0.0")
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

    fout = open(args.out_filename,"wb")

    if args.feature_matrix != None:
        feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)
        fout.write("id1\tscore\t%s\n" % "\t".join(args.header_names))
    else:
        fout.write("id1\tscore\n")

    for clustid, cluster in enumerate(clusters):
        for prot_pair in it.combinations(cluster,2):
            id1 = "%s_%s" % (clustid, prot_pair[0])
            id2 = "%s_%s" % (clustid, prot_pair[1])
            try:
                score = pairwise_interactions[frozenset([prot_pair[0],prot_pair[1]])]
            except KeyError:
                continue
            if args.feature_matrix != None:
                out_list = []
                ids_str_order = "%s" % sorted([prot_pair[0],prot_pair[1]])
                for field in args.features:
                    out_list.append(str(feature_table[feature_table['frozenset_ids_str_order'] == ids_str_order][field].values[0] != 0.0))

                fout.write("%s (pp) %s\t%s\t%s" % (id1, id2, score, "\t".join(out_list)))
            else:
                fout.write("%s (pp) %s\t%s" % (id1, id2, score))
            fout.write("\n")

    fout.close()

if __name__ == "__main__":
    main()
