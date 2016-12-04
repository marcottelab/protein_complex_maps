
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Outputs table of protein complex edges from database")
    parser.add_argument("--outputfile", action="store", dest="outputfile", required=False, default=None,
                                    help="Filename to store output, default=stdout")
    parser.add_argument("--field_delimiter", action="store", dest="field_delimiter", required=False, default='\t',
                                    help="delimiter between fields, default = \\t")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    complexes = db.session.query(cdb.Complex).all()

    output_str = "complex_id%scorr_pval%st_count%sq_count%sqandt_count%sqandt_by_q%sqandt_by_t%sterm_id%st_type%st_group%st_name%sdepth_in_group%sqandt_list\n" % tuple([args.field_delimiter]*12)

    for comp in complexes:
        for enrichment in comp.enrichments:

            output_str = output_str + "%s"*25 % (
                                                        comp.complex_id, args.field_delimiter, 
                                                        enrichment.corr_pval, args.field_delimiter,
                                                        enrichment.t_count, args.field_delimiter,  
                                                        enrichment.q_count, args.field_delimiter,  
                                                        enrichment.qandt_count, args.field_delimiter,  
                                                        enrichment.qandt_by_q, args.field_delimiter,  
                                                        enrichment.qandt_by_t, args.field_delimiter,  
                                                        enrichment.term_id, args.field_delimiter,  
                                                        enrichment.t_type, args.field_delimiter,  
                                                        enrichment.t_group, args.field_delimiter,  
                                                        enrichment.t_name, args.field_delimiter,  
                                                        enrichment.depth_in_group, args.field_delimiter,  
                                                        enrichment.qandt_list 
                                                        ) + "\n"
    if args.outputfile != None:
        f = open(args.outputfile, "wb")
        f.write(output_str)
        f.close()
    else:
        print output_str

if __name__ == "__main__":
    main()


