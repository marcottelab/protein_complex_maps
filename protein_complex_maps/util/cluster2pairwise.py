

import argparse
import itertools as it

def main():

    parser = argparse.ArgumentParser(description="Tool to convert clusters to pairs with option of adding cluster id to output pairs ")
    parser.add_argument("--filename", action="store", dest="filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
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
            if args.add_cluster_id:
                if set(["%s_%s" % (clustid,prot_pair[0]),"%s_%s" % (clustid,prot_pair[1])]) not in output_sets:
                    fout.write("%s_%s" % (clustid, prot_pair[0]))
                    fout.write("\t")
                    fout.write("%s_%s" % (clustid, prot_pair[1]))
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




