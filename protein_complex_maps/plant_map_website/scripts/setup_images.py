from __future__ import print_function
import argparse
import numpy as np
import itertools as it

import csv

import sys


#import protein_complex_maps.plant_map_website.complex_db as cdb


# Make consistent with other loader scripts
def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--cluster_table", action="store", dest="node_table", required=True, 
                                    help="Filename node table")

    args = parser.parse_args()


    node_table_file = open(args.node_table,"rb")
    csvreader = csv.reader(node_table_file, delimiter=',', quotechar='"')
    count = 1
 
    with open("clustid_key.csv", "w") as outfile:
       outfile.write("OrthogroupID,clustid,clustid_set,order\n")
       for line in csvreader:

        # CDM header
        # ID, cut_0.99, cut_0.95, cut_0.9, cut_0.85, cut_0.8, cut_0.75, cut_0.7, cut_0.65, cut_0.6, cut_0.55, cut_0.5, cut_0.45, cut_0.4, cut_0.35, cut_0.3, cut_0.25, cut_0.2, cut_0.15, cut_0.1
        # Select cut_0.7  cut_0.6  cut_0.45 cut_0.25cut_0.7  cut_0.6  cut_0.45 cut_0.25

        #kdrew: header: acc,clustid,clustid_key,genename,key,proteinname,uniprot_link
        #kdrew: if header do not parse
            if 'dendrogram_order' in line:
                continue
        #print line8,10,13, 17
            OrthogroupID = line[0]
            clustid_1 = line[7]
            clustid_2 = line[9]
            clustid_3 = line[12]
            clustid_4 = line[16]
            order = line[-1].strip()

            outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_1, "clustid_1", order))
            outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_2, "clustid_2", order))
            outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_3, "clustid_3", order))
            outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_4, "clustid_4", order))

        #print("prior", OrthogroupID, clustid_1, clustid_2, clustid_3, clustid_4)

if __name__ == "__main__":
    main()


