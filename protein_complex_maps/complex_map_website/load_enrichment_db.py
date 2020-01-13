
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads enrichment sql tables from input files")
    parser.add_argument("--enrichment_file", action="store", dest="enrichment_file", required=True, 
                                    help="Filename enrichment table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    enrichment_df = pd.read_csv(args.enrichment_file)
    for _, row in enrichment_df.iterrows():

        #index,complex_id,description,effective_domain_size,goshv,group_id,intersection_genes,intersection_size,name,native,p_value,parents,precision,query,query_size,recall,significant,source,source_order,term_size

        c = db.session.query(cdb.Complex).filter_by(complex_id=row.complex_id).first()
        if c:
            print "complex id: %s" % c.id
            ce = cdb.get_or_create(db, cdb.ComplexEnrichment, 
                complex_key           = c.id,
                description           = row.description,
                effective_domain_size = row.effective_domain_size,
                goshv                 = row.goshv,                
                group_id              = row.group_id,             
                intersection_genes    = row.intersection_genes,    
                intersection_size     = row.intersection_size,    
                #name                  = row.name, #kdrew: name is overloaded here and defaults to name of row rather than the column 'name'                
                name                  = row['name'],                 
                native                = row.native,               
                p_value               = row.p_value,              
                parents               = row.parents,              
                precision             = row.precision,            
                query                 = row.query,                
                query_size            = row.query_size,           
                recall                = row.recall,               
                significant           = row.significant,          
                source                = row.source,               
                source_order          = row.source_order,         
                term_size             = row.term_size)             

            db.session.add(ce)
            db.session.commit()
        else:
            print "Cannot find complex %s" % (complex_id)


if __name__ == "__main__":
    main()


