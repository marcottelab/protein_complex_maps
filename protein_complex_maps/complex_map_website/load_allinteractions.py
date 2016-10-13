
import argparse
import numpy as np
import itertools as it

import csv

#import protein_complex_maps.complex_map_website.complex_db as cdb
import complex_db as cdb


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
    for line in interaction_table_file.readlines():

        #if 'complex:' in line#:
        #    complex_id = line.split()[1]
        #    continue

        #kdrew: if header do not parse
        #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list
        #if 'signf' in line:
        #    continue

        print line
        split_line = line.split('\t')
        print split_line
        prot1 = split_line[0]
        prot2 = split_line[1]
        svmscore = float(split_line[2].strip())
        p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot1).first()
        p2 = db.session.query(cdb.Protein).filter_by(gene_id=prot2).first()
        if p1 and p2:
            print "protein id1: %s" % p1.id
            print "protein id2: %s" % p2.id
            edge = cdb.get_or_create(db, cdb.Edge,
                                        protein_key = p1.id,
                                        protein_key2 = p2.id,
                                        score = score,
                                        )



            db.session.add(edge)
            db.session.commit()
        else:
            print "Cannot find complex %s" % (complex_id)


if __name__ == "__main__":
    main()


