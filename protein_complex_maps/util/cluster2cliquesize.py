

import argparse
import itertools as it

def main():

    parser = argparse.ArgumentParser(description="Tool to convert clusters to cliques")
    parser.add_argument("--filename", action="store", dest="filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
                                            help="Output filename ")
    parser.add_argument("--clique_size", action="store", type=int, dest="clique_size", required=False, default=2,
                                            help="Size of clique, default=2 (pairwise)")

    args = parser.parse_args()

    #kdrew: error checking
    if not args.filename:
        print "\nError: Specify --filename \n"
        parser.print_help()
        return

    clusters = []
    f = open(args.filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())

    #kdrew: dumps all pairs in clusters one per line in entrez id
    fout = open(args.out_filename,"wb")

    #kdrew: used to store pairs that have already been outputted, removes redundancy
    output_sets = []
    for clustid, cluster in enumerate(clusters):
        for cliqueid, prot_set in enumerate(it.combinations(cluster,args.clique_size)):
            for prot_pair in it.combinations(prot_set,2):
                if set(["%s_%s_%s" % (clustid,cliqueid, prot_pair[0]),"%s_%s_%s" % (clustid,cliqueid, prot_pair[1])]) not in output_sets:
                    fout.write("%s_%s_%s" % (clustid, cliqueid, prot_pair[0]))
                    fout.write("\t")
                    fout.write("%s_%s_%s" % (clustid, cliqueid, prot_pair[1]))
                    fout.write("\n")
                    output_sets.append(set(["%s_%s_%s" % (clustid, cliqueid, prot_pair[0]),"%s_%s_%s" % (clustid, cliqueid, prot_pair[1])]))
    fout.close()

if __name__ == "__main__":
    main()




