
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Outputs table of protein complexes from database")
    parser.add_argument("--outputfile", action="store", dest="outputfile", required=False, default=None,
                                    help="Filename to store output, default=stdout")
    parser.add_argument("--field_delimiter", action="store", dest="field_delimiter", required=False, default=',',
                                    help="delimiter between fields, default = ,")
    parser.add_argument("--infield_delimiter", action="store", dest="infield_delimiter", required=False, default='\t',
                                    help="delimiter within fields, default = \\t")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    complexes = db.session.query(cdb.Complex).all()

    output_str = "complex_id%sgene_ids%sgenenames\n" % (args.field_delimiter,args.field_delimiter)
    for comp in complexes:
        print comp.complex_id
        geneid_str = args.infield_delimiter.join([p.gene_id for p in comp.proteins])
        genename_str = args.infield_delimiter.join([p.genename() for p in comp.proteins])
        output_str = output_str + "%s%s%s%s%s\n" % (comp.complex_id, args.field_delimiter,geneid_str,args.field_delimiter,genename_str)

    if args.outputfile != None:
        f = open(args.outputfile, "wb")
        f.write(output_str)
        f.close()
    else:
        print output_str


if __name__ == "__main__":
    main()


