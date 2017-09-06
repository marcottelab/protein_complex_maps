
import argparse
import numpy as np
import itertools as it

import csv
import matplotlib.pyplot as plt
#import protein_complex_maps.complex_map_website.complex_db as cdb
import complex_db as cdb
import pandas as pd

def make_data_frame(query, columns):
    """
    Takes a sqlalchemy query and a list of columns, returns a dataframe.
    """
    def make_row(x):
        return dict([(c, getattr(x, c)) for c in columns])       
    return pd.DataFrame([make_row(x) for x in query])

def main():

    parser = argparse.ArgumentParser(description="Loads all edge interactions sql tables from input files")

    parser.add_argument("--all_interactions", action="store", dest="interaction_file", required=True,
                                    help="Filename pairwise interaction table with scores")



    #cdm: db_data/blake_bioplex_prey_hein_revisitTrain_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt



    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    interaction_table_file = open(args.interaction_file,"rb")
    complex_id = None
    count = 1
    for line in interaction_table_file.readlines():
        if count % 100 == 0:
            print count
        count = count + 1    
        #if 'complex:' i        #print line
        split_line = line.split(' ')
        #print split_line
        group1 = split_line[0]
        group2 = split_line[1]
        score = float(split_line[2].strip())
        edge = cdb.get_or_create(db, cdb.Edge,
                                        GroupID_key = group1,
                                        GroupID_key2 = group2,
                                        in_complex = 0, 
                                        score = score,
                                        )
        db.session.add(edge)
        db.session.commit()

if __name__ == "__main__":
    main()


