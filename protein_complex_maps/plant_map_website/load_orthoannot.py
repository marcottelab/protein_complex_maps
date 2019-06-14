
import argparse
import numpy as np
import itertools as it

import csv
import matplotlib.pyplot as plt
#import protein_complex_maps.complex_map_website.complex_db as cdb
import complex_db as cdb
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Loads a annotation table")

    parser.add_argument("--annotation_file", action="store", dest="annotation_file", required=True,
                                    help="At least 3 column tab-separated file of eggnogID Arath_genename, orthoannotation")


    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    annotation_table = open(args.annotation_file,"rb")
    complex_id = None
    count = 1

#ID      arath_genenames arath_Entries   arath_Entry_names       arath_Protein_names     disruptions     tair_disruptionslloyd2012_LOFs  arath_functions arath_misc      pathway unipathway      BioCyc  Reactome        BRENDA  kegg_pws       ec       arath_masses    arath_protein_names     arath_GO        devstages       tissues tair    araport Annotation     orysj_genenames  orysj_Entries   orysj_Entry_names       orysj_Protein_names     orysj_disruptions       orysj_functions
#orysj_misc


    with(open(args.annotation_file,"rb")) as annotation_table:
        
        for line in annotation_table.readlines():
            if count % 100 == 0:
                print count
            count = count + 1    
            #if 'complex:' i        #print line
            if "arath_genenames" in line:
                print(line)
                continue 
            split_line = line.split('\t')
            #print(split_line)
            try:
                OrthogroupID = split_line[0]
                ArathGenenames = split_line[1]
                EggnogAnnot = split_line[24]
                Tair = split_line[23].strip("\n")
            #print(OrthogroupID, ArathGenenames, EggnogAnnot, Tair)
                a = cdb.get_or_create(db, cdb.Orthoannot, EggnogAnnot = EggnogAnnot, Tair = Tair, ArathGenenames = ArathGenenames)
                o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
                db.session.add(a)
                db.session.commit()
    
                oam = cdb.get_or_create(db, cdb.OrthogroupAnnotMapping, orthogroup_key=o.id, orthoannot_key=a.id)
                db.session.add(oam)
                db.session.commit()    
            except Exception as e:
                 continue

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

