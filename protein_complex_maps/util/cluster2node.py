import argparse
import itertools as it
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")

    args = parser.parse_args()

    protid_set = set()

    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())
        map( protid_set.add, line.split() )

    d = dict()
    d['clustid'] = []
    d['ID'] = []
    d['clustid_key'] = []

    for clustid, cluster in enumerate(clusters):
        for group_id in cluster:
            clust_id_key = "%s_%s" % (clustid, group_id)

            d['clustid'].append(clustid)
            d['ID'].append(str(group_id))
            d['clustid_key'].append(clust_id_key)

    df = pd.DataFrame(d)
    df.to_csv(args.output_filename)


if __name__ == "__main__":
    main()
