

import argparse
import itertools as it

import protein_complex_maps.protein_util as pu

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

    input_ids = list(set([x for cluster in clusters for x in cluster]))
    inputID2ACC_map = pu.map_protein_ids(input_ids, args.from_id, "ACC", reviewed=args.reviewed)
    flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    ACC2outputID_map = pu.map_protein_ids(flatten_list, "ACC", args.to_id, reviewed=args.reviewed)

    if not args.pairwise:
        #kdrew: dumps clusters one per line in output id
        fout = open(args.out_filename,"wb")
        for cluster in clusters:
            cluster_prot = []
            for prot in cluster:
                #print prot
                out_id = None
                for acc in inputID2ACC_map[prot]:
                    if len(ACC2outputID_map[acc]) == 0:
                        continue
                    else:
                        out_id = ACC2outputID_map[acc]
                        break
                if  out_id != None and len(out_id) > 0:
                    cluster_prot.append(out_id[0])

                    #fout.write(out_id[0])
                    #fout.write("\t")
                
            fout.write("\t".join(cluster_prot))
            fout.write("\n")
        fout.close()


    else:
        #kdrew: dumps all pairs in clusters one per line in entrez id
        fout = open(args.out_filename,"wb")

        #kdrew: used to store pairs that have already been outputted, removes redundancy
        output_sets = []
        for cluster in clusters:
            for prot_pair in it.combinations(cluster,2):
                #print prot
                out_id1 = None
                out_id2 = None
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
                    fout.write(out_id1[0])
                    fout.write("\t")
                    fout.write(out_id2[0])
                    fout.write("\n")
                    output_sets.append(set([out_id1[0],out_id2[0]]))
        fout.close()

if __name__ == "__main__":
    main()




