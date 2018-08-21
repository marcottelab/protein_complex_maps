
import argparse
import numpy as np
import itertools as it

import csv
import pandas as pd

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--node_table", action="store", dest="node_table", required=True, 
                                    help="Filename node table")
    parser.add_argument("--humap_corum_mapping", action="store", dest="humap_corum_mapping", required=True,
                                    help="Filename HuMap Corum mapping")

    #kdrew: /home/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_prey_hein_prey_revisitTrain/trim_subunits/blake_bioplex_prey_hein_prey_revisitTrain_corum_train_allComplexesCore_trainSplit_noTestOverlap_psweep7.ii149.clusterone_agglomod.ii94.reduced.trimThreshold_reduced.nodeTable.txt

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()
    
    df = pd.read_csv(args.node_table)
    df_mapping = pd.read_csv(args.humap_corum_mapping)

    print df.head()

    for index, row in df.iterrows():
        p = cdb.get_or_create(db, cdb.Protein, uniprot_acc=row.ID, genename=row['Gene names  (primary )'], proteinname=row['Protein names'])
        db.session.add(p)

        corum_id = df_mapping.query("humap_corum_id == @row.clustid")['corum_ids'].values[0]
        humap_id = df_mapping.query("humap_corum_id == @row.clustid")['humap_ids'].values[0]
        print corum_id
        print humap_id                                                            

        rnp_label = row.RNP_label
        if rnp_label != rnp_label:
            rnp_label = ''

        c = cdb.get_or_create(db, cdb.Complex, complex_id=row.clustid, rnp_label=rnp_label, corum_id=corum_id, humap_id=humap_id)
        db.session.add(c)
        db.session.commit()

        print "complex id: %s" % c.id
        print "protein id: %s" % p.id
        pcm = cdb.get_or_create(db, cdb.ProteinComplexMapping, protein_key=p.id, complex_key=c.id)
        db.session.add(pcm)
        db.session.commit()


if __name__ == "__main__":
    main()


