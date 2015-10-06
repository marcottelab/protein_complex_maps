

import argparse
import itertools as it

import protein_complex_maps.protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table to pairwise interactions (and scores) with corresponding cluster_ids" )
    parser.add_argument("--pairwise_filename", action="store", dest="pairwise_filename", required=True,
                                            help="Filename of pairwise file, (ie. one line per ppi, third column contains score)")
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
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

    fout.write("id1\tscore\n")
    for clustid, cluster in enumerate(clusters):
        for prot_pair in it.combinations(cluster,2):
            id1 = "%s_%s" % (clustid, prot_pair[0])
            id2 = "%s_%s" % (clustid, prot_pair[1])
            try:
                score = pairwise_interactions[frozenset([prot_pair[0],prot_pair[1]])]
            except KeyError:
                continue
            fout.write("%s (pp) %s\t%s" % (id1, id2, score))
            fout.write("\n")

    fout.close()

if __name__ == "__main__":
    main()
