
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    #kdrew: should really generalize the bait prey files and make them optional, also need to reorganize the readin through pandas
    #kdrew: but just trying to get the db loaded

    parser = argparse.ArgumentParser(description="Loads edge sql tables from input files")
    parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
                                    help="Filename edge table")
    parser.add_argument("--bioplex_bait_prey_file", action="store", dest="bioplex_bait_prey_file", required=True, 
                                    help="Filename of bioplex bait prey pairs (format: bait_geneid, prey_geneid)")
    parser.add_argument("--hein_bait_prey_file", action="store", dest="hein_bait_prey_file", required=True, 
                                    help="Filename of hein bait prey pairs (format: bait_geneid, prey_geneid)")
    parser.add_argument("--boldt_bait_prey_file", action="store", dest="boldt_bait_prey_file", required=True, 
                                    help="Filename of boldt bait prey pairs (format: bait_geneid, prey_geneid)")
    parser.add_argument("--youn_bait_prey_file", action="store", dest="youn_bait_prey_file", required=True, 
                                    help="Filename of youn bait prey pairs (format: bait_geneid, prey_geneid)")
    parser.add_argument("--gupta_bait_prey_file", action="store", dest="gupta_bait_prey_file", required=True, 
                                    help="Filename of gupta bait prey pairs (format: bait_geneid, prey_geneid)")

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

    boldt_bait_preys = set()
    boldt_bait_prey_file = open(args.boldt_bait_prey_file,"rb")
    for line in boldt_bait_prey_file.readlines():
        boldt_bait_preys.add(tuple([x.strip() for x in line.split(',')]))

    youn_bait_preys = set()
    youn_bait_prey_file = open(args.youn_bait_prey_file,"rb")
    for line in youn_bait_prey_file.readlines():
        youn_bait_preys.add(tuple([x.strip() for x in line.split(',')]))

    gupta_bait_preys = set()
    gupta_bait_prey_file = open(args.gupta_bait_prey_file,"rb")
    for line in gupta_bait_prey_file.readlines():
        gupta_bait_preys.add(tuple([x.strip() for x in line.split(',')]))

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
        id1 = split_line[1]
        prot1 = id1.split('(pp)')[0].split('_')[1].strip()
        prot2 = id1.split('(pp)')[1].split('_')[1].strip()
        score = float(split_line[2])
        evidence_dict = dict()
        #id1     score   fractions       bioplex hein    bioplex_prey    hein_prey       Guru    Malo    bioplex2        gupta_ciliated  gupta_nonciliated       boldt   youn    bioplex2_hygeo  gupta_hygeo     boldt_hygeo     youn_hygeo      treiber_hygeo   hygeo_only
        evidence_dict['fraction'] = ('True' in split_line[3])
        evidence_dict['bioplex'] = ('True' in split_line[4])
        evidence_dict['hein'] = ('True' in split_line[5])
        evidence_dict['bioplex_WMM'] = ('True' in split_line[6])
        evidence_dict['hein_WMM'] = ('True' in split_line[7])
        evidence_dict['Guru'] = ('True' in split_line[8])
        evidence_dict['Malo'] = ('True' in split_line[9])
        evidence_dict['bioplex2'] = ('True' in split_line[10])
        evidence_dict['gupta_ciliated'] = ('True' in split_line[11])
        evidence_dict['gupta_nonciliated'] = ('True' in split_line[12])
        evidence_dict['boldt'] = ('True' in split_line[13])
        evidence_dict['youn'] = ('True' in split_line[14])
        evidence_dict['bioplex2_WMM'] = ('True' in split_line[15])
        evidence_dict['gupta_WMM'] = ('True' in split_line[16])
        evidence_dict['boldt_WMM'] = ('True' in split_line[17])
        evidence_dict['youn_WMM'] = ('True' in split_line[18])
        evidence_dict['treiber_WMM'] = ('True' in split_line[19])
        evidence_dict['WMM_only'] = ('True' in split_line[20])
    

        p1 = db.session.query(cdb.Protein).filter_by(gene_id=prot1).first()
        p2 = db.session.query(cdb.Protein).filter_by(gene_id=prot2).first()

        if evidence_dict['bioplex'] or evidence_dict['bioplex2']:
            bioplex_baits = []
            if tuple([prot1,prot2]) in bioplex_bait_preys:
                bioplex_baits.append(p1.genename())
            if tuple([prot2,prot1]) in bioplex_bait_preys:
                bioplex_baits.append(p2.genename())
            bioplex_baits.sort()
            #print bioplex_baits

        if evidence_dict['hein']:
            hein_baits = []
            if tuple([prot1,prot2]) in hein_bait_preys:
                hein_baits.append(p1.genename())
            if tuple([prot2,prot1]) in hein_bait_preys:
                hein_baits.append(p2.genename())
            hein_baits.sort()

        if evidence_dict['boldt']:
            boldt_baits = []
            if tuple([prot1,prot2]) in boldt_bait_preys:
                boldt_baits.append(p1.genename())
            if tuple([prot2,prot1]) in boldt_bait_preys:
                boldt_baits.append(p2.genename())
            boldt_baits.sort()

        if evidence_dict['youn']:
            youn_baits = []
            if tuple([prot1,prot2]) in youn_bait_preys:
                youn_baits.append(p1.genename())
            if tuple([prot2,prot1]) in youn_bait_preys:
                youn_baits.append(p2.genename())
            youn_baits.sort()

        if evidence_dict['gupta_nonciliated'] or evidence_dict['gupta_ciliated']:
            gupta_baits = []
            if tuple([prot1,prot2]) in gupta_bait_preys:
                gupta_baits.append(p1.genename())
            if tuple([prot2,prot1]) in gupta_bait_preys:
                gupta_baits.append(p2.genename())
            gupta_baits.sort()

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
                    if k == 'bioplex' or k == 'bioplex2':
                        kstr = 'bioplex (%s)' % (','.join(bioplex_baits))
                    if k == 'bioplex_WMM' or k == 'bioplex2_WMM':
                        kstr = 'bioplex_WMM' 
                    if k == 'boldt':
                        if len(boldt_baits) > 0:
                            kstr+= ' (%s)' % (','.join(boldt_baits))
                    if k == 'youn':
                        kstr+= ' (%s)' % (','.join(youn_baits))
                    if k == 'gupta_ciliated' or k == 'gupta_nonciliated':
                        kstr = 'gupta (%s)' % (','.join(gupta_baits))

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


