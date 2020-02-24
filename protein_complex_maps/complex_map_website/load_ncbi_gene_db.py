
import argparse
import numpy as np
import pandas as pd
import itertools as it

import csv

import protein_complex_maps.util.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads gene sql tables ")
    parser.add_argument("--ncbi_file", action="store", dest="ncbi_file", required=True, 
                            help="Filename ncbi table (gene2annotation file from ftp://ftp.ncbi.nih.gov/gene/DATA/)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    ncbi_df = pd.read_table(args.ncbi_file)


    #kdrew: get all proteins in database
    proteins = db.session.query(cdb.Protein).filter_by(uniprot_acc = '').all()
    #proteins = db.session.query(cdb.Protein).filter_by(gene_id=80199).all()
    geneid_set = set([p.gene_id for p in proteins])
    
    #for prot_id in protid_set:
    #    print prot_id
    #    p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot_id).first()
    #    print "id: %s" % p1.id

    #kdrew: retreive from uniprot all accs
    #inputID2ACC_map = pu.map_protein_ids(list(protid_set), "P_ENTREZGENEID", "ACC", reviewed=True)
    #flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    #kdrew: retreive from uniprot all genenames
    #kdrew: added return_list=True so all genenames are returned rather than just the first
    #genename_map = pu.get_from_uniprot(flatten_list, 'genes', return_list=True)
    #print genename_map

    #kdrew: for each geneid in database
    for gene_id1 in geneid_set:

        #kdrew: get protein with that geneid
        p1 = db.session.query(cdb.Protein).filter_by(gene_id=gene_id1).first()
        if p1:
            #kdrew: add all genenames to Gene table and link to protein instance
            for genename in ncbi_df.query("GeneID == @gene_id1")['Symbol'].unique():
                print "%s %s" % (genename, gene_id1)
                if genename != None:
                    #gene = cdb.get_or_create(db,cdb.Gene, gene_id = gene_id1, genename=genename, protein_key = p1.id)
                    gene = cdb.get_or_create(db,cdb.Gene, genename=genename, protein_key = p1.id)
                    db.session.add(gene)
                    db.session.commit()
        else:
            print "Cannot find proteins %s" % (gene_id1)


if __name__ == "__main__":
    main()


