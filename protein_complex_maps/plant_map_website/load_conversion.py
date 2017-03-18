
import argparse
import numpy as np
import itertools as it

import csv
import matplotlib.pyplot as plt
#import protein_complex_maps.complex_map_website.complex_db as cdb
import complex_db as cdb
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Loads a conversion table")

    parser.add_argument("--conversion_file", action="store", dest="conversion_file", required=True,
                                    help="Four column tab-separated table of genename, proteinname, gene_id, uniprot_acc")



    #cdm: db_data/blake_bioplex_prey_hein_revisitTrain_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt



    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    conversion_table = open(args.conversion_file,"rb")
    complex_id = None
    count = 1
    for line in conversion_table.readlines():
        if count % 100 == 0:
            print count
        count = count + 1    
        #if 'complex:' i        #print line
        split_line = line.split(',')
        genename = split_line[0]
        proteinname = split_line[1]
        gene_id = split_line[2]
        uniprot_acc = split_line[3].strip()
        convert = cdb.get_or_create(db, cdb.Conversion,
                                        genename = genename,
                                        proteinname = proteinname,
                                        gene_id = gene_id,
                                        uniprot_acc = uniprot_acc,
                                        )
        db.session.add(convert)
        db.session.commit()
if __name__ == "__main__":

    main()


