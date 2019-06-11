
import argparse
import numpy as np
import itertools as it

import csv

import sys

print(sys.path)

import protein_complex_maps.plant_map_website.complex_db as cdb

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
    for line in csvreader:
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

        print("prior", OrthogroupID, clustid_1, clustid_2, clustid_3, clustid_4)

        h = cdb.get_or_create(db, cdb.Hiercomplex, clustid_1 = clustid_1, clustid_2 = clustid_2, clustid_3 = clustid_3, clustid_4 = clustid_4)
        print("what")
        o = cdb.get_or_create(db, cdb.Orthogroup)
        db.session.add(h)
        db.session.add(o)
        db.session.commit()

        print "complex id: %s" % h.id
        print "orthogroup id: %s" % o.id
        ocm = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h.id)
        db.session.add(ocm)
        db.session.commit()


if __name__ == "__main__":
    main()


