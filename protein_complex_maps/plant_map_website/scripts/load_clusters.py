from __future__ import print_function
import argparse
import numpy as np
import itertools as it

import csv

import sys


#import protein_complex_maps.plant_map_website.complex_db as cdb
import complex_db as cdb


# Make consistent with other loader scripts
def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--cluster_table", action="store", dest="node_table", required=True, 
                                    help="Filename node table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    node_table_file = open(args.node_table,"rb")
    csvreader = csv.reader(node_table_file, delimiter=',', quotechar='"')
    count = 1
    for line in csvreader:
        if count % 100 == 0:
            print(count)
        count = count + 1

        # CDM header
        # ID, cut_0.99, cut_0.95, cut_0.9, cut_0.85, cut_0.8, cut_0.75, cut_0.7, cut_0.65, cut_0.6, cut_0.55, cut_0.5, cut_0.45, cut_0.4, cut_0.35, cut_0.3, cut_0.25, cut_0.2, cut_0.15, cut_0.1
        # Select cut_0.7  cut_0.6  cut_0.45 cut_0.25cut_0.7  cut_0.6  cut_0.45 cut_0.25

        #kdrew: header: acc,clustid,clustid_key,genename,key,proteinname,uniprot_link
        #kdrew: if header do not parse
        if 'ID' in line:
            continue
        #print line8,10,13, 17
        OrthogroupID = line[0]
        clustid_1 = line[7]
        clustid_2 = line[9]
        clustid_3 = line[12]
        clustid_4 = line[16].strip()

        #print("prior", OrthogroupID, clustid_1, clustid_2, clustid_3, clustid_4)

        h1 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_1, clustid_set = 'clustid_1') 
        h2 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_2, clustid_set = 'clustid_2') 
        h3 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_3, clustid_set = 'clustid_3') 
        h4 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_4, clustid_set = 'clustid_4') 

        o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
        db.session.add(h1)
        db.session.add(h2)
        db.session.add(h3)
        db.session.add(h4)

        db.session.add(o)
        db.session.commit()

        ocm1 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h1.id)
        ocm2 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h2.id)
        ocm3 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h3.id)
        ocm4 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h4.id)
 

        db.session.add(ocm1)
        db.session.add(ocm2)
        db.session.add(ocm3)
        db.session.add(ocm4)

        db.session.commit()


if __name__ == "__main__":
    main()


