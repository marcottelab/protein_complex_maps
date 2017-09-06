
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
    parser.add_argument("--field_delimiter", action="store", dest="field_delimiter", required=False, default=',',
                                    help="delimiter between fields, default = ,")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    edges = db.session.query(cdb.Edge).all()

    proteins = set()
    protein_edges = set()
    output_str = "gene_id1%sgene_id2%sgenename1%sgenename2%sscore\n" % (args.field_delimiter,args.field_delimiter,args.field_delimiter,args.field_delimiter)

    for edge in edges:
        prots = edge.get_proteins()
        output_str = output_str + "%s%s%s%s%s%s%s%s%s\n" % (prots[0].gene_id, args.field_delimiter, 
                                                            prots[1].gene_id, args.field_delimiter,
                                                            prots[0].genename(), args.field_delimiter,
                                                            prots[1].genename(), args.field_delimiter,
                                                            edge.score)
        proteins.add(prots[0].gene_id)
        proteins.add(prots[1].gene_id)
        protein_edges.add(frozenset([prots[0].gene_id,prots[1].gene_id]))
        
    if args.outputfile != None:
        f = open(args.outputfile, "wb")
        f.write(output_str)
        f.close()
    else:
        print output_str

    print "#ofProteins: %s, #ofEdges: %s" % (len(proteins),len(protein_edges))

if __name__ == "__main__":
    main()


