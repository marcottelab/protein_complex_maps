#import sys
import argparse
import itertools as it
import pandas as pd
#sys.path.append('/project/cmcwhite/protein_complex_maps/protein_complex_maps')
#import protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")
    parser.add_argument("--annotation_filename", action="store", dest="annotation_filename", required=True,
                                            help="Annotation. 1 row per ID, IDs match network")

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
    print(df)
    annot=pd.read_csv(args.annotation_filename, index_col=False, sep="\t")

    annot = annot.set_index(["ID"])
    print(annot)
    print(df)
    final = df.join(annot, how="left", on="ID")

    final.to_csv(args.output_filename)


if __name__ == "__main__":
    main()
