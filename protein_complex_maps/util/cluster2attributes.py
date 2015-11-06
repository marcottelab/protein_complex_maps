

import argparse
import itertools as it

import protein_complex_maps.protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to get genenames" )
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

    attributes = dict()
    attributes['genes'] = dict()
    #attributes['protein names'] = dict()
    #kdrew: used to store pairs that have already been outputted, removes redundancy
    input_ids = list(set([x for cluster in clusters for x in cluster]))
    #results = pu.get_from_uniprot(input_ids, attributes.keys(), return_as_string=True)
    results = pu.get_from_uniprot(input_ids, 'genes')
    print results

    for clustid, cluster in enumerate(clusters):
        for protid in cluster:
            if args.add_cluster_id:
                attributes['genes']["%s_%s" % (clustid, protid)] = results[protid]
            else:
                attributes['genes'][protid] = results[protid]

    #kdrew: write header
    fout.write("protid\t%s\n" % ('\t'.join(attributes.keys())))
    for protid in attributes['genes']:
        fout.write("%s\t%s" % (protid, attributes['genes'][protid]))
        fout.write("\n")

        #for prot_pair in it.combinations(cluster,2):
        #    if args.add_cluster_id:
        #        if set(["%s_%s" % (clustid,prot_pair[0]),"%s_%s" % (clustid,prot_pair[1])]) not in output_sets:
        #            fout.write("%s_%s" % (clustid, prot_pair[0]))
        #            fout.write("\t")
        #            fout.write("%s_%s" % (clustid, prot_pair[1]))
        #            fout.write("\n")
        #            output_sets.append(set(["%s_%s" % (clustid, prot_pair[0]),"%s_%s" % (clustid, prot_pair[1])]))
        #    else:
        #        if set([prot_pair[0],prot_pair[1]]) not in output_sets:
        #            fout.write(prot_pair[0])
        #            fout.write("\t")
        #            fout.write(prot_pair[1])
        #            fout.write("\n")
        #            output_sets.append(set([prot_pair[0],prot_pair[1]]))
    fout.close()

if __name__ == "__main__":
    main()
