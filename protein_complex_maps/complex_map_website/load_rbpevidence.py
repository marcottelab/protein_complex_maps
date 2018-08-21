
import argparse
import numpy as np
import itertools as it

import csv
import pandas as pd

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--rbpevidence_table", action="store", dest="rbpevidence_table", required=True, 
                                    help="Filename of RBP evidence table for individual proteins")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()
    
    df = pd.read_csv(args.rbpevidence_table).set_index('ACC')
    print df.head()
    print df.dtypes

    for index, row in df.iterrows():
        #print row
        p = db.session.query(cdb.Protein).filter_by(uniprot_acc=index).first()


        if p == None:
            continue

        #ACC,low_throughput,high_throughput,computational,rnabinding,rnabinding_uniprot,trendel,queiroz,huang,bao,HEK293-RIC_Hs_Baltz2012,HuH7-IC_Hs_Beckmann2015,HeLa-RNPxl_Hs_Kramer2014,HeLa-IC_Hs_Castello2012,HeLa-RBDmap_Hs_Castello2016,K562-serIC-chr_Hs_Conrad2016_chr,K562-serIC_Hs_Conrad2016,Enzyme,Metabolism,"""Metabolic.Enzyme""",annotation

        if row['trendel']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Trendel et al. 2018')
            db.session.add(rbpe)
        if row['queiroz']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Queiroz et al. 2018')
            db.session.add(rbpe)
        if row['huang']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Huang et al. 2018')
            db.session.add(rbpe)
        if row['bao']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Bao et al. 2018')
            db.session.add(rbpe)
        if row['HEK293-RIC_Hs_Baltz2012']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Baltz et al. 2012')
            db.session.add(rbpe)
        if row['HuH7-IC_Hs_Beckmann2015']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Beckmann et al. 2015')
            db.session.add(rbpe)
        if row['HeLa-RNPxl_Hs_Kramer2014']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Kramer et al. 2014')
            db.session.add(rbpe)
        if row['HeLa-IC_Hs_Castello2012']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Castello et al. 2012')
            db.session.add(rbpe)
        if row['HeLa-RBDmap_Hs_Castello2016']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Castello et al. 2016')
            db.session.add(rbpe)
        if row['K562-serIC-chr_Hs_Conrad2016_chr']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Conrad et al. 2016_chr')
            db.session.add(rbpe)
        if row['K562-serIC_Hs_Conrad2016']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='high throughput', evidence='Conrad et al. 2016')
            db.session.add(rbpe)

        #if row['Keywords'] == row['Keywords'] and 'Ribonucleoprotein' in row['Keywords']:
        #    rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='Uniprot RNP')
        #    db.session.add(rbpe)
        #if row['rnabinding']:
        #    rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='low throughput', evidence_type='GO RNA Binding')
        #    db.session.add(rbpe)
        if row['rnabinding']:
            rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='low throughput', evidence='GO RNA Binding')
            db.session.add(rbpe)

        #if row['sliding_pvalues_fdrcor'] <= 0.05:
        #    rbpe = cdb.get_or_create(db, cdb.RBPEvidence, protein_key=p.id, evidence_type='diffrac', evidence='DIFFRAC (significant)')
        #    db.session.add(rbpe)

        db.session.commit()

if __name__ == "__main__":
    main()


