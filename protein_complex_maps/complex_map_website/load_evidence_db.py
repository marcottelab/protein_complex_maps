
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads edge sql tables from input files")
    parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
                                    help="Filename edge table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    edge_table_file = open(args.edge_file,"rb")
    for line in edge_table_file.readlines():

        #kdrew: if header do not parse
        #id1     score   fractions       bioplex hein    bioplex_prey    hein_prey
        if 'score' in line:
            continue

        print line
        split_line = line.split('\t')
        print split_line

        #kdrew: id1 example: Q96BY6  Q9Y3L3  0.999999
        prot1 = split_line[0]
        prot2 = split_line[1]
        score = float(split_line[2])

        p1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=prot1).first()
        p2 = db.session.query(cdb.Protein).filter_by(uniprot_acc=prot2).first()

        if p1 and p2:
            #kdrew: enforce order on protein ids
            if p2.id < p1.id:
                p2, p1 = p1, p2
            print "protein id1: %s" % p1.id
            print "protein id2: %s" % p2.id
            edge = cdb.get_or_create(db, cdb.Edge, 
                                        protein_key = p1.id,
                                        protein_key2 = p2.id,
                                        score = score,
                                        )
            db.session.add(edge)
            db.session.commit()

        else:
            print "Cannot find proteins %s" % (id1)


if __name__ == "__main__":
    main()


