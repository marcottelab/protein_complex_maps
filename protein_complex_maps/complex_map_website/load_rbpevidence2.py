
import argparse
import numpy as np
import itertools as it

import csv
import pandas as pd

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--rbpstats_table", action="store", dest="rbpstats_table", required=True, 
                                    help="Filename of RBP stats table for individual proteins")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()
    
    df = pd.read_csv(args.rbpstats_table, index_col=0)
    print df.head()

    for index, row in df.iterrows():
        #print row
        p = db.session.query(cdb.Protein).filter_by(uniprot_acc=index).first()


        if p == None:
            continue

        if row['Keywords'] == row['Keywords'] and 'Ribonucleoprotein' in row['Keywords']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='low throughput', evidence='Uniprot RNP')
            db.session.add(rbpe)

        if row['sliding_pvalues_fdrcor'] <= 0.05:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='diffrac', evidence='DIFFRAC (significant)')
            db.session.add(rbpe)

        db.session.commit()

if __name__ == "__main__":
    main()


