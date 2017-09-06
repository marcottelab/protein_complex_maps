
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads edge sql tables from input files")
    parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
                                    help="Filename edge table")
    parser.add_argument("--bioplex_bait_prey_file", action="store", dest="bioplex_bait_prey_file", required=True, 
                                    help="Filename of bioplex bait prey pairs (format: bait_geneid, prey_geneid)")
    parser.add_argument("--hein_bait_prey_file", action="store", dest="hein_bait_prey_file", required=True, 
                                    help="Filename of hein bait prey pairs (format: bait_geneid, prey_geneid)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    bioplex_bait_preys = set()
    bioplex_bait_prey_file = open(args.bioplex_bait_prey_file,"rb")
    for line in bioplex_bait_prey_file.readlines():
        bioplex_bait_preys.add(tuple([x.strip() for x in line.split(',')]))

    hein_bait_preys = set()
    hein_bait_prey_file = open(args.hein_bait_prey_file,"rb")
    for line in hein_bait_prey_file.readlines():
        hein_bait_preys.add(tuple([x.strip() for x in line.split(',')]))

    print bioplex_bait_preys
    print hein_bait_preys

    edge_table_file = open(args.edge_file,"rb")
    for line in edge_table_file.readlines():

        #kdrew: if header do not parse
        #id1     score   fractions       bioplex hein    bioplex_prey    hein_prey
        if 'score' in line:
            continue

        print line
        split_line = line.split('\t')
        print split_line

        #kdrew: id1 example: 0_153129 (pp) 0_10670 
        id1 = split_line[0]
        prot1 = id1.split('(pp)')[0].split('_')[1].strip()
        prot2 = id1.split('(pp)')[1].split('_')[1].strip()
        score = float(split_line[1])
        evidence_dict = dict()
        evidence_dict['fraction'] = ('True' in split_line[2])
        evidence_dict['bioplex'] = ('True' in split_line[3])
        evidence_dict['hein'] = ('True' in split_line[4])
        evidence_dict['bioplex_prey'] = ('True' in split_line[5])
        evidence_dict['hein_prey'] = ('True' in split_line[6])

        p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot1).first()
        p2 = db.session.query(cdb.Protein).filter_by(gene_id=prot2).first()

        if evidence_dict['bioplex']:
            bioplex_baits = []
            if tuple([prot1,prot2]) in bioplex_bait_preys:
                bioplex_baits.append(p1.genename())
            if tuple([prot2,prot1]) in bioplex_bait_preys:
                bioplex_baits.append(p2.genename())
            bioplex_baits.sort()

        if evidence_dict['hein']:
            hein_baits = []
            if tuple([prot1,prot2]) in hein_bait_preys:
                hein_baits.append(p1.genename())
            if tuple([prot2,prot1]) in hein_bait_preys:
                hein_baits.append(p2.genename())
            hein_baits.sort()

        if p1 and p2:
            #kdrew: enforce order on protein ids
            if p2.id < p1.id:
                p2, p1 = p1, p2
            print "protein id1: %s" % p1.id
            print "protein id2: %s" % p2.id
            edge = cdb.get_or_create(db, cdb.Edge, 
                                        protein_key = p1.id,
                                        protein_key2 = p2.id,
                                        score = score,
                                        )
            db.session.add(edge)
            db.session.commit()

            for k in evidence_dict:
                if evidence_dict[k]:
                    kstr = k
                    if k == 'hein':
                        kstr+= ' (%s)' % (','.join(hein_baits))
                    if k == 'bioplex':
                        kstr+= ' (%s)' % (','.join(bioplex_baits))
                    evidence = cdb.get_or_create(db, cdb.Evidence,
                                                    edge_key = edge.id,
                                                    evidence_type = kstr
                                                    )
                    db.session.add(evidence)
                    db.session.commit()

        else:
            print "Cannot find proteins %s" % (id1)


if __name__ == "__main__":
    main()


