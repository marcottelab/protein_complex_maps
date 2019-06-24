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

    parser = argparse.ArgumentParser(description="Loads a conversion table")

    parser.add_argument("--conversion_file", action="store", dest="conversion_file", required=True,
                                    help="At least 3 column comma-separated file of eggnogID Spec ProteinID")



    #cdm: db_data/blake_bioplex_prey_hein_revisitTrain_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt



    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    conversion_table = open(args.conversion_file,"rb")
    complex_id = None
    count = 1
    t0 = time.time()
    opms = []
    with(open(args.conversion_file,"rb")) as conversion_table:

        for line in conversion_table.readlines():
            if "ProteinID" in line:
                 continue
            if count % 1000 == 0:
                print(count, str(time.time() - t0))
                t0 = time.time()
            count = count + 1
            #if 'complex:' i        #print line
            split_line = line.split(',')
            OrthogroupID = split_line[0]
            ProteinID = split_line[1]
            Spec = split_line[2].strip("\n")
            #convert = cdb.get_or_create(db, cdb.Conversion,
            #                                OrthogroupID = OrthogroupID,
            #                                Spec = Spec,
            #                                ProteinID = ProteinID,
            #                                )

            o = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID).first()
            p = db.session.query(cdb.Protein).filter_by(ProteinID = ProteinID).first()
            #db.session.flush()
            opm = cdb.get_or_create(db, cdb.OrthogroupProteinMapping, orthogroup_key=o.id, protein_key=p.id)
            opms.append(opm)
      
    db.session.add_all(opms)
    db.session.commit()


if __name__ == "__main__":
    main()
