
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--node_table", action="store", dest="node_table", required=True, 
                                    help="Filename node table")

    #kdrew: /home/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_prey_hein_prey_revisitTrain/trim_subunits/blake_bioplex_prey_hein_prey_revisitTrain_corum_train_allComplexesCore_trainSplit_noTestOverlap_psweep7.ii149.clusterone_agglomod.ii94.reduced.trimThreshold_reduced.nodeTable.txt

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    node_table_file = open(args.node_table,"rb")
    csvreader = csv.reader(node_table_file, delimiter=',', quotechar='"')
    for line in csvreader:
        #kdrew: header: acc,clustid,clustid_key,genename,key,proteinname,uniprot_link
        #kdrew: if header do not parse
        if 'clustid_key' in line:
            continue
        #print line
        acc = line[0]
        clustid = line[1]
        clustid_key = line[2]
        genename = line[3]
        geneid = line[4]
        proteinname = line[5].strip()
        uniprot_link = line[6]

        p = cdb.get_or_create(db, cdb.Protein, gene_id = geneid, uniprot_acc=acc, genename=genename, proteinname=proteinname, uniprot_url=uniprot_link)
        db.session.add(p)

        if clustid != '':
            c = cdb.get_or_create(db, cdb.Complex, complex_id=clustid)
            db.session.add(c)
            db.session.commit()

            print "complex id: %s" % c.id
            print "protein id: %s" % p.id
            pcm = cdb.get_or_create(db, cdb.ProteinComplexMapping, protein_key=p.id, complex_key=c.id)
            db.session.add(pcm)
            db.session.commit()


if __name__ == "__main__":
    main()


