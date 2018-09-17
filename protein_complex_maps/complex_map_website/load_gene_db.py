
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.util.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads gene sql tables ")
    #parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
    #                                help="Filename edge table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()


    #kdrew: get all proteins in database
    proteins = db.session.query(cdb.Protein).all()
    print len(proteins)
    #proteins = db.session.query(cdb.Protein).filter_by(gene_id=80199).all()
    protid_set = set([p.uniprot_acc for p in proteins])
    
    #for prot_id in protid_set:
    #    print prot_id
    #    p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot_id).first()
    #    print "id: %s" % p1.id

    #kdrew: retreive from uniprot all accs
    #inputID2ACC_map = pu.map_protein_ids(list(protid_set), "P_ENTREZGENEID", "ACC", reviewed=True)
    #flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    #kdrew: retreive from uniprot all genenames
    #kdrew: added return_list=True so all genenames are returned rather than just the first
    genename_map = pu.get_from_uniprot(list(protid_set), 'genes', return_list=True)
    print genename_map

    #kdrew: for each acc in database
    for prot_id in protid_set:

        #kdrew: get protein with that acc
        p1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=prot_id).first()
        if p1:
            #kdrew: add all genenames to Gene table and link to protein instance
            try:
                for genename in genename_map[prot_id]:
                    print genename
                    if genename != None:
                        gene = cdb.get_or_create(db,cdb.Gene, genename=genename, protein_key = p1.id)
                        db.session.add(gene)
                        db.session.commit()
            except KeyError:
                print "No genename for %s" % (prot_id)
        else:
            print "Cannot find proteins %s" % (prot_id)


if __name__ == "__main__":
    main()


