

import argparse
import itertools as it

def main():

    parser = argparse.ArgumentParser(description="Tool to convert clusters to pairs with option of adding cluster id to output pairs ")
    parser.add_argument("--filename", action="store", dest="filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--remove_overlap_filename", action="store", dest="remove_overlap_filename", required=False, default=None,
                                            help="Filename of cluster file, (ie. one line per cluster), edges in this set will not be included in the final set")
    parser.add_argument("--keep_overlap_filename", action="store", dest="keep_overlap_filename", required=False, default=None,
                                            help="Filename of cluster file, (ie. one line per cluster), edges in this set will only be included in the final set")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
                                            help="Output filename ")
    parser.add_argument("--add_cluster_id", action="store_true", dest="add_cluster_id", required=False, default=False,
                                            help="Add a cluster id to the protein id. Used to distinguish same protein in multiple clusters, default=False")

    args = parser.parse_args()

    #kdrew: error checking
    if not args.filename:
        print "\nError: Specify --filename \n"
        parser.print_help()
        return

    #kdrew: generate set of pairs in the removal set if filename is given on commandline
    removal_pairs_set = set()
    if args.remove_overlap_filename != None:
        #kdrew: read in removal_pairs
        removal_f = open(args.remove_overlap_filename, "rb")
        for line in removal_f.readlines():
            for pair in it.combinations( set( line.split() ), 2):
                removal_pairs_set.add(frozenset(pair))

    #kdrew: generate set of pairs in the keep set if filename is given on commandline
    keep_pairs_set = set()
    #kdrew: initialize keep file as the input file in case there is no keep overlap file (defaults to keep all pairs)
    ko_filename = args.filename
    if args.keep_overlap_filename != None:
        ko_filename = args.keep_overlap_filename

    #kdrew: read in keep_pairs
    keep_f = open(ko_filename, "rb")
    for line in keep_f.readlines():
        for pair in it.combinations( set( line.split() ), 2):
            keep_pairs_set.add(frozenset(pair))


    clusters = []
    f = open(args.filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())

    #kdrew: dumps all pairs in clusters one per line in entrez id
    fout = open(args.out_filename,"wb")

    #kdrew: used to store pairs that have already been outputted, removes redundancy
    output_sets = []
    for clustid, cluster in enumerate(clusters):
        for prot_pair in it.combinations(cluster,2):
            if frozenset(prot_pair) not in removal_pairs_set and frozenset(prot_pair) in keep_pairs_set:
                if args.add_cluster_id:
                    if set(["%s_%s" % (clustid,prot_pair[0]),"%s_%s" % (clustid,prot_pair[1])]) not in output_sets:
                        fout.write("%s_%s" % (clustid, prot_pair[0]))
                        fout.write("\t")
                        fout.write("%s_%s" % (clustid, prot_pair[1]))
                        fout.write("\t")
                        fout.write("%s\t%s" % (prot_pair[0], prot_pair[1]))
                        fout.write("\n")
                        output_sets.append(set(["%s_%s" % (clustid, prot_pair[0]),"%s_%s" % (clustid, prot_pair[1])]))
                else:
                    if set([prot_pair[0],prot_pair[1]]) not in output_sets:
                        fout.write(prot_pair[0])
                        fout.write("\t")
                        fout.write(prot_pair[1])
                        fout.write("\n")
                        output_sets.append(set([prot_pair[0],prot_pair[1]]))
    fout.close()

if __name__ == "__main__":
    main()




