from __future__ import print_function
import argparse
#import numpy as np

#import itertools as it

#import csv
import time
#import numpy as np
#import itertools as it

#import csv
#import matplotlib.pyplot as plt
#import protein_complex_maps.complex_map_website.complex_db as cdb
import complex_db as cdb
#import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Loads a score table")

    parser.add_argument("--score_file", action="store", dest="score_file", required=True,
                                    help="3 column tab-separated file of eggnogID eggnogID score")



    #cdm: db_data/blake_bioplex_prey_hein_revisitTrain_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt



    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    score_table = open(args.score_file,"rb")
    count = 1
    t0 = time.time()
    scrs = []
    with(open(args.score_file,"rb")) as score_table:

        for line in score_table.readlines():
            if "P_1" in line: #header
                 continue
            if count % 1000 == 0:
                print(count, str(time.time() - t0))
                t0 = time.time()

            count = count + 1

            split_line = line.split(',')
            OrthogroupID1 = split_line[0]
            OrthogroupID2 = split_line[1]
            ScoreVal = split_line[2].strip("\n")


            o1 = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID1).first()
            o2 = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID2).first()
            #db.session.flush()
            scr = cdb.get_or_create(db, cdb.Score, OrthogroupID1_key=o1.id, OrthogroupID2_key=o2.id, ScoreVal = ScoreVal)
            scrs.append(scr)
      
    db.session.add_all(scrs)
    db.session.commit()


if __name__ == "__main__":
    main()
