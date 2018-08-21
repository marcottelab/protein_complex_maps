
import argparse
import numpy as np
import itertools as it

import csv
import pandas as pd

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--rnpstats_table", action="store", dest="rnpstats_table", required=True, 
                                    help="Filename of RNP stats table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()
    
    df = pd.read_csv(args.rnpstats_table)
    print df.head()

    for index, row in df.iterrows():
        print row
        c = db.session.query(cdb.Complex).filter_by(complex_id=index).first()

        if c == None:
            continue

        crs = cdb.get_or_create(db, cdb.ComplexRNPStats, complex_key=c.id, 
                                                        sec2_cntl_median_correlation         =row['sec2_cntl_median_correlation'],
                                                        sec2_rnase_median_correlation        =row['sec2_rnase_median_correlation'],
                                                        sec2_cntl_num_ids                    =row['sec2_cntl_num_ids'],
                                                        sec2_cntl_fraction_ids               =row['sec2_cntl_fraction_ids'],
                                                        sec2_rnase_num_ids                   =row['sec2_rnase_num_ids'],
                                                        sec2_rnase_fraction_ids              =row['sec2_rnase_fraction_ids'],
                                                        sec2_num_fdr05                       =row['sec2_num_fdr05'],
                                                        sec2_fraction_fdr05                  =row['sec2_fraction_fdr05'],
                                                        sec2_num_fdr2                        =row['sec2_num_fdr2'],
                                                        sec2_fraction_fdr2                   =row['sec2_fraction_fdr2'],
                                                        sec2_num_annotated                   =row['sec2_num_annotated'],
                                                        sec2_fraction_annotated              =row['sec2_fraction_annotated'],
                                                        sec2_num_rnabinding                  =row['sec2_num_rnabinding'],
                                                        sec2_fraction_rnabinding             =row['sec2_fraction_rnabinding'],
                                                        sec2_num_rnabinding_uniprot          =row['sec2_num_rnabinding_uniprot'],
                                                        sec2_fraction_rnabinding_uniprot     =row['sec2_fraction_rnabinding_uniprot'],
                                                        sec2_num_lt                          =row['sec2_num_lt'],
                                                        sec2_fraction_lt                     =row['sec2_fraction_lt'],
                                                        sec2_num_ht                          =row['sec2_num_ht'],
                                                        sec2_fraction_ht                     =row['sec2_fraction_ht'],
                                                        sec2_num_unannotated                 =row['sec2_num_unannotated'],
                                                        sec2_fraction_unannotated            =row['sec2_fraction_unannotated'],
                                                        sec2_num_ht_go_uniprot_sec2_rbp      =row['sec2_num_ht_go_uniprot_sec2_rbp'],
                                                        sec2_fraction_ht_go_uniprot_sec2_rbp =row['sec2_fraction_ht_go_uniprot_sec2_rbp'],
                                                        sec2_logratio_psm                    =row['sec2_logratio_psm'],
                                                        sec2_bw_exps_median_correlation      =row['sec2_bw_exps_median_correlation'],
                                                        sec2_num_fdr5                        =row['sec2_num_fdr5'],
                                                        sec2_fraction_fdr5                   =row['sec2_fraction_fdr5'],
                                                        Significant_RNP                      =row['Significant_RNP'],
                                                        Significant_20FDR_08Corr_RNP         =row['Significant_20FDR_08Corr_RNP'],
                                                        High_Throughput_RNP                  =row['High_Throughput_RNP'],
                                                        GO_Uniprot_RNP                       =row['GO_Uniprot_RNP'] )

        db.session.add(crs)
        db.session.commit()


if __name__ == "__main__":
    main()


