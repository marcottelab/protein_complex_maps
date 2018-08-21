
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

        rbps = cdb.get_or_create(db, cdb.RBPStats, protein_key=p.id, 

        #,diffrac,diffrac_percent,diffrac_normalized,emd,pearsonr,mean_abundance,annotated,zscore,sliding_zscore,pvalues,pvalues_fdrcor,sliding_pvalues,sliding_pvalues_fdrcor,low_throughput,high_throughput,computational,rnabinding,rnabinding_uniprot,trendel,queiroz,huang,bao,HEK293-RIC_Hs_Baltz2012,HuH7-IC_Hs_Beckmann2015,HeLa-RNPxl_Hs_Kramer2014,HeLa-IC_Hs_Castello2012,HeLa-RBDmap_Hs_Castello2016,K562-serIC-chr_Hs_Conrad2016_chr,K562-serIC_Hs_Conrad2016,Enzyme,Metabolism,"""Metabolic.Enzyme""",annotation,Entry name,Status,Protein names,Gene names,Organism,Cross-reference (GeneID),Gene names  (primary ),Keywords,Gene ontology (biological process),Gene ontology (cellular component),Gene ontology (molecular function),ht_go_uniprot_sec2_rbp

        diffrac                   = row['diffrac'],                 
        diffrac_percent           = row['diffrac_percent'],          
        diffrac_normalized        = row['diffrac_normalized'],       
        emd                       = row['emd'],                     
        pearsonr                  = row['pearsonr'],             
        mean_abundance            = row['mean_abundance'],           
        zscore                    = row['zscore'],                   
        sliding_zscore            = row['sliding_zscore'],           
        pvalues                   = row['pvalues'],                  
        pvalues_fdrcor            = row['pvalues_fdrcor'],           
        sliding_pvalues           = row['sliding_pvalues'],          
        sliding_pvalues_fdrcor    = row['sliding_pvalues_fdrcor'])

        db.session.add(rbps)
        db.session.commit()

if __name__ == "__main__":
    main()


