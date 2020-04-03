
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.util.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads uniprot entries into database")
    parser.add_argument("--uniprot_file", action="store", dest="uniprot_file", required=True, 
                            help="Uniprot filename (format: ACC Protein_name    GeneID  primary_genename    genenames)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    f = open(args.uniprot_file,"rb")
    for line in f.readlines():
        lsplit = line.split("\t")
        ACC = lsplit[0]
        pname = lsplit[1]
        geneids = lsplit[2].split(';')
        genenames = lsplit[4].split()
        uniprot_url = "http://www.uniprot.org/uniprot/%s" % (ACC)

        for geneid in geneids:
            if geneid != "":
                protein = cdb.get_or_create(db,cdb.Protein, gene_id = geneid, uniprot_acc = ACC, proteinname = pname, uniprot_url = uniprot_url)
                db.session.add(protein)
                db.session.commit()

                for genename in genenames:
                    #kdrew: do a little clean up
                    gname = genename.split(';')[0]
                    gene = cdb.get_or_create(db,cdb.Gene, genename=gname, protein_key = protein.id)
                    db.session.add(gene)
                    db.session.commit()


        #kdrew: some entries in uniprot annotation file do not have geneid so cannot map protein in database exactly, just search based on uniprot acc instead
        if len(geneids) == 1:
            protein = db.session.query(cdb.Protein).filter_by(uniprot_acc=ACC).first()
            if protein:
                for genename in genenames:
                    #kdrew: do a little clean up
                    gname = genename.split(';')[0]
                    gene = cdb.get_or_create(db,cdb.Gene, genename=gname, protein_key = protein.id)
                    db.session.add(gene)
                    db.session.commit()



    #kdrew: get all proteins in database
    #proteins = db.session.query(cdb.Protein).all()
    #protid_set = set([p.gene_id for p in proteins])
    

#    #kdrew: for each geneid in database
#    for prot_id in protid_set:
#
#        #kdrew: get protein with that geneid
#        p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot_id).first()
#        if p1:
#            #kdrew: add all genenames to Gene table and link to protein instance
#            for acc in inputID2ACC_map[prot_id]:
#                for genename in genename_map[acc]:
#                    print genename
#                    if genename != None:
#                        gene = cdb.get_or_create(db,cdb.Gene, gene_id = prot_id, genename=genename, protein_key = p1.id)
#                        db.session.add(gene)
#                        db.session.commit()
#        else:
#            print "Cannot find proteins %s" % (prot_id)
#

if __name__ == "__main__":
    main()


