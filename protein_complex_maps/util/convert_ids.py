

import argparse
import itertools as it

import protein_complex_maps.util.protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to convert ids from one type to another for cluster files")
    parser.add_argument("--filename", action="store", dest="filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
                                            help="Output filename ")
    parser.add_argument("--from_id", action="store", dest="from_id", required=True, 
                                            help="input id type, (ENSEMBL_ID,etc)")
    parser.add_argument("--to_id", action="store", dest="to_id", required=False, default="ACC", 
                                            help="convert ids to this type, (P_ENTREZGENEID,ACC,etc), default=ACC")
    parser.add_argument("--reviewed", action="store_true", dest="reviewed", required=False, default=False,
                                            help="map only to reviewed ids, default=False")
    parser.add_argument("--pairwise", action="store_true", dest="pairwise", required=False, default=False,
                                            help="Split complex up into all pairs")
    parser.add_argument("--add_cluster_id", action="store_true", dest="add_cluster_id", required=False, default=False,
                                            help="Add a cluster id to the protein id. Used to distinguish same protein in multiple clusters, default=False")
    parser.add_argument("--ppi_with_scores", action="store_true", dest="ppi_with_scores", required=False, default=False,
                                            help="Input filename is pairwise with score in 3rd column, default=False")
    parser.add_argument("--columns", action="store", dest="columns", nargs='+', required=False, default=[],
                                            help="Convert ids in specified columns, outputting the remaining columns unchanged")
    parser.add_argument("--orig_id", action="store_true", dest="orig_id", required=False, default=False,
                                            help="Adds the original id instead of None if conversion fails")
    parser.add_argument("--contact_email", action="store", dest="contact_email", required=True, 
                                            help="This script queries uniprot webservice which requires a contact email")

    args = parser.parse_args()

    #kdrew: error checking
    if not args.filename:
        print "\nError: Specify --filename \n"
        parser.print_help()
        return

    clusters = []
    #kdrew: save so you can write out ppi score later
    ppi_scores = dict()
    additional_columns = dict()
    f = open(args.filename, "rb")
    for i, line in enumerate(f.readlines()):
        if args.ppi_with_scores:
            clusters.append(line.split()[:2])
            ppi_scores[i] = line.split()[2]
        #kdrew: untested
        elif len(args.columns) > 0:
            clust = []
            additional_c = []
            for c in range(len(line.split())):
                if str(c) in args.columns:
                    clust.append(line.split()[c])
                else:
                    additional_c.append(line.split()[c])
            clusters.append(clust)
            additional_columns[i] = additional_c
        else:
            clusters.append(line.split())

    input_ids = list(set([x for cluster in clusters for x in cluster]))
    inputID2ACC_map = pu.map_protein_ids(input_ids, args.from_id, "ACC", args.contact_email, reviewed=args.reviewed)
    flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    ACC2outputID_map = pu.map_protein_ids(flatten_list, "ACC", args.to_id, args.contact_email, reviewed=args.reviewed)

    if not args.pairwise:
        output_string = ""
        #kdrew: dumps clusters one per line in output id
        for clustid, cluster in enumerate(clusters):
            cluster_prot = []
            for prot in cluster:
                #print prot

                out_id = None
                if args.orig_id:
                    #kdrew: default id is the original id
                    out_id = [prot,]

                for acc in inputID2ACC_map[prot]:
                    if len(ACC2outputID_map[acc]) == 0:
                        continue
                    else:
                        out_id = ACC2outputID_map[acc]
                        break
                if  out_id != None and len(out_id) > 0:
                    if args.add_cluster_id:
                        cluster_prot.append("%s_%s" % (clustid, out_id[0]))
                    else:
                        cluster_prot.append(out_id[0])

                    #fout.write(out_id[0])
                    #fout.write("\t")
                
            #fout.write("\t".join(cluster_prot))
            output_string += "\t".join(cluster_prot)
            if args.ppi_with_scores:
                #fout.write("\t%s" % ppi_scores[clustid])
                output_string += "\t%s" % ppi_scores[clustid]
            #fout.write("\n")
            elif len(args.columns) > 0:
                output_string += "\t%s" % ("\t".join(additional_columns[clustid]))


            output_string += "\n"

        fout = open(args.out_filename,"wb")
        fout.write(output_string)
        fout.close()


    else:
        output_string = ""
        #kdrew: used to store pairs that have already been outputted, removes redundancy
        output_sets = []
        for clustid, cluster in enumerate(clusters):
            for prot_pair in it.combinations(cluster,2):
                #print prot

                out_id1 = None
                out_id2 = None
                if args.orig_id:
                    #kdrew: default id is the original id
                    out_id1 = [prot_pair[0],]
                    out_id2 = [prot_pair[1],]

                for acc1 in inputID2ACC_map[prot_pair[0]]:
                    if len(ACC2outputID_map[acc1]) == 0:
                        continue
                    else:
                        out_id1 = ACC2outputID_map[acc1]
                        break
                for acc2 in inputID2ACC_map[prot_pair[1]]:
                    if len(ACC2outputID_map[acc2]) == 0:
                        continue
                    else:
                        out_id2 = ACC2outputID_map[acc2]
                        break

                if out_id1 != None and out_id2 != None and len(out_id1) > 0 and len(out_id2) > 0 and set([out_id1[0],out_id2[0]]) not in output_sets:
                    if args.add_cluster_id:
                        #fout.write("%s_%s" % (clustid, out_id1[0]))
                        #fout.write("\t")
                        #fout.write("%s_%s" % (clustid, out_id2[0]))

                        output_string += "%s_%s" % (clustid, out_id1[0])
                        output_string += "\t"
                        output_string += "%s_%s" % (clustid, out_id2[0])
                        if args.ppi_with_scores:
                            #fout.write("\t%s" % ppi_scores[clustid])
                            output_string += "\t%s" % ppi_scores[clustid]
                        elif len(args.columns) > 0:
                            output_string += "\t%s" % ("\t".join(additional_columns[clustid]))
                            
                        output_string += "\n"

                        output_sets.append(set(["%s_%s" % (clustid, out_id1[0]),"%s_%s" % (clustid, out_id2[0])]))
                    else:
                        #fout.write(out_id1[0])
                        #fout.write("\t")
                        #fout.write(out_id2[0])
                        output_string += out_id1[0]
                        output_string += "\t"
                        output_string += out_id2[0]
                        if args.ppi_with_scores:
                            #fout.write("\t%s" % ppi_scores[clustid])
                            output_string += "\t%s" % ppi_scores[clustid]
                        elif len(args.columns) > 0:
                            output_string += "\t%s" % ("\t".join(additional_columns[clustid]))
                        #fout.write("\n")
                        output_string += "\n"
                        output_sets.append(set([out_id1[0],out_id2[0]]))

        #kdrew: dumps all pairs in clusters one per line in entrez id
        fout = open(args.out_filename,"wb")
        fout.write(output_string)
        fout.close()

if __name__ == "__main__":
    main()




