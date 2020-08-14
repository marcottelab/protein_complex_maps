
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Outputs complex database to flat files")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    fout = open(args.output_file,"w")
    fout.write("HuMAP2_ID,Confidence,Uniprot_ACCs,genenames\n")
    for c in db.session.query(cdb.Complex).all():
        if c:
            print "complex id: %s" % c.humap2_id
            print "proteins: %s" % ' '.join([p.uniprot_acc for p in c.proteins])
            print "genenames: %s" % ' '.join([p.genename() for p in c.proteins])
            if len(c.proteins) == 0:
                continue
            fout.write("%s,%s,%s,%s\n" % (c.humap2_id, c.top_rank, ' '.join([p.uniprot_acc for p in c.proteins]), ' '.join([p.genename() for p in c.proteins])))
        else:
            print "Cannot find complex %s" % (complex_id)

    fout.close()

if __name__ == "__main__":
    main()


