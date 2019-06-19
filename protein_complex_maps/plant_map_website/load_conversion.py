from __future__ import print_function
import argparse
#import numpy as np

#import itertools as it

import csv
import time
#from time import sleep
import complex_db as cdb
import pandas as pd


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
    t0  = time.time()
    with(open(args.conversion_file,"rb")) as conversion_table:
       
        os = []
        ps = [] 
        opms = [] 
        for line in conversion_table.readlines():
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
            # Find a better way to do this. Avoid concurrent loading and error
            error_count = 0
            try:    
                o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
                p = cdb.get_or_create(db, cdb.Protein, ProteinID = ProteinID, Spec = Spec)
                #db.session.add(o)
                #db.session.add(p)
                #db.session.commit()
        
                #opm = cdb.get_or_create(db, cdb.OrthogroupProteinMapping, orthogroup_key=o.id, protein_key=p.id)
                #db.session.add(opm)
                #db.session.commit()
            except Exception as E:
                error_count = error_count + 1
                if error_count < 1:
                    sleep(1)                   
                    o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
                    p = cdb.get_or_create(db, cdb.Protein, ProteinID = ProteinID, Spec = Spec)
                    os.append(o)
                    ps.append(p)
                    #db.session.add(o)
                    #db.session.add(p)
            
                    #opm = cdb.get_or_create(db, cdb.OrthogroupProteinMapping, orthogroup_key=o.id, protein_key=p.id)
                    #db.session.add(opm)
                    #opms.append(opm)
                    #db.session.commit()
                else: 
                   print(E)
                   continue
               
            if count % 100 == 0:
                print(count, str(time.time() - t0))
                t0 = time.time()
                db.session.add_all(os)
                db.session.add_all(ps)
                for i in range(len(os)):
                   opm = cdb.get_or_create(db, cdb.OrthogroupProteinMapping, orthogroup_key=os[i].id, protein_key=ps[i].id)
                   opms.append(opm)
                db.session.add_all(opms)
                os = []
                ps = []
                opms = []
                db.session.flush()


    db.session.commit()

  

if __name__ == "__main__":

    main()
#with open('million_users.csv', 'r') as csv_file:
#    csv_reader = csv.reader(csv_file)
#
#buffer = []
#for row in reader:
#    buffer.append({
#        'OrthogroupI': row[0],
#         'ProteinID': row[1],
#          'Spec': row[2]
#    })
#    if len(buffer) % 10000 == 0:
#        session.bulk_insert_mappings(buffer)
#        buffer = []
#
#session.bulk_insert_mappings(buffer)

