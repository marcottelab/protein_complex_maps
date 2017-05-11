
import argparse
import numpy as np
import itertools as it

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

    enrichment_table_file = open(args.enrichment_file,"rb")
    complex_id = None
    for line in enrichment_table_file.readlines():

        if 'complex:' in line:
            complex_id = line.split()[1]
            continue

        #kdrew: if header do not parse
        #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list
        if 'signf' in line:
            continue

        print line
        split_line = line.split('\t')
        print split_line

        corr_pval = split_line[2]
        t_count = split_line[3] 
        q_count = split_line[4] 
        qandt_count = split_line[5]  
        qandt_by_q = split_line[6] 
        qandt_by_t = split_line[7] 
        term_id = split_line[8] 
        t_type = split_line[9] 
        t_group = split_line[10] 
        t_name = split_line[11] 
        depth_in_group = split_line[12] 
        qandt_list = split_line[13].strip()

        c = db.session.query(cdb.Complex).filter_by(complex_id=complex_id).first()
        if c:
            print "complex id: %s" % c.id
            ce = cdb.get_or_create(db, cdb.ComplexEnrichment, 
                                        complex_key=c.id, 
                                        corr_pval=corr_pval, 
                                        t_count=t_count, 
                                        q_count=q_count, 
                                        qandt_count=qandt_count, 
                                        qandt_by_q=qandt_by_q, 
                                        qandt_by_t=qandt_by_t, 
                                        term_id=term_id, 
                                        t_type=t_type, 
                                        t_group=t_group, 
                                        t_name=t_name, 
                                        depth_in_group=depth_in_group, 
                                        qandt_list=qandt_list)
            db.session.add(ce)
            db.session.commit()
        else:
            print "Cannot find complex %s" % (complex_id)


if __name__ == "__main__":
    main()


