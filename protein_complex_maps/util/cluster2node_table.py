

import argparse
import itertools as it
import pandas as pd

import protein_complex_maps.protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")
    parser.add_argument("--from_id", action="store", dest="from_id", required=True, 
                                            help="input id type, (P_ENTREZGENEID, ENSEMBL_ID,etc)")
    parser.add_argument("--reviewed", action="store_true", dest="reviewed", required=False, default=False,
                                            help="map only to reviewed ids, default=False")
    args = parser.parse_args()

    protid_set = set()

    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())
        map( protid_set.add, line.split() )

    inputID2ACC_map = pu.map_protein_ids(list(protid_set), args.from_id, "ACC", reviewed=args.reviewed)
    flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    genename_map = pu.get_from_uniprot(flatten_list, 'genes')

    d = dict()
    d['clustid'] = []
    d['key'] = []
    d['clustid_key'] = []
    d['acc'] = []
    d['genename'] = []

    for clustid, cluster in enumerate(clusters):
        for prot_id in cluster:
            clust_id_key = "%s_%s" % (clustid, prot_id)

            d['clustid'].append(clustid)
            d['key'].append(prot_id)
            d['clustid_key'].append(clust_id_key)
            try:
                d['acc'].append(inputID2ACC_map[prot_id][0])
                d['genename'].append(genename_map[inputID2ACC_map[prot_id][0]])
            except IndexError:
                d['genename'].append(None)



    df = pd.DataFrame(d)
    df.to_csv(args.output_filename, index=False)


if __name__ == "__main__":
    main()
