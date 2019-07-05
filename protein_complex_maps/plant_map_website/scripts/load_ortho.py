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

    parser = argparse.ArgumentParser(description="Loads a orthogroup table")

    parser.add_argument("--orthogroup_file", action="store", dest="orthogroup_file", required=True,
                                    help="One column of _unique_ orthogroupIDs")



    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    orthogroup_table = open(args.orthogroup_file,"rb")
    count = 1
    t0  = time.time()
    with(open(args.orthogroup_file,"rb")) as orthogroup_table:
       
        os = []




        for line in orthogroup_table.readlines():
            count = count + 1    
            OrthogroupID = line.strip("\n")
            #o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
            o = cdb.Orthogroup(OrthogroupID = OrthogroupID) # Plain load
            os.append(o)  
            if count % 10000 == 0:
                 
                print(count, str(time.time() - t0))
                t0 = time.time()
                #print(os, ps)
    db.session.add_all(os)
    db.session.flush()  # Gets each object updated with its .id
    print("added all IDs!")
    db.session.commit()

  

if __name__ == "__main__":

    main()
