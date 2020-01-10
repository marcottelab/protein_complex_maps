
import argparse
import numpy as np
import itertools as it
from sqlalchemy import func, or_, and_

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Populates EdgeComplexMapping table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    #kdrew: iterate through complexes
    complexes = db.session.query(cdb.Complex).all()
    for c in complexes:
        for prot1, prot2 in it.combinations(c.proteins,2):
            #kdrew: edge table enforces order
            if prot2.id < prot1.id:
                prot2, prot1 = prot1, prot2

            e = db.session.query(cdb.Edge).filter( and_(cdb.Edge.protein_key == prot1.id, cdb.Edge.protein_key2 == prot2.id) ).first()
            if e != None:
                print "complex id: %s" % c.id
                print "edge id: %s" % e.id
                ecm = cdb.get_or_create(db, cdb.EdgeComplexMapping, edge_key=e.id, complex_key=c.id)
                db.session.add(ecm)
                db.session.commit()


if __name__ == "__main__":
    main()


