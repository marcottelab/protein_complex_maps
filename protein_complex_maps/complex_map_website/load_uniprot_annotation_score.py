
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads annotation score from input files")
    parser.add_argument("--annotation_score_file", action="store", dest="annotation_score_file", required=True, 
                                    help="Uniprot annotation score filename")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    score_df = pd.read_table(args.annotation_score_file)
    for _, row in score_df.iterrows():

        p = db.session.query(cdb.Protein).filter_by(uniprot_acc=row.Entry).first()
        if p:
            print "complex id: %s" % p.id
            p.annotation_score = row.Annotation
            db.session.commit()
        else:
            print "Cannot find protein %s" % (row.Entry)


if __name__ == "__main__":
    main()


