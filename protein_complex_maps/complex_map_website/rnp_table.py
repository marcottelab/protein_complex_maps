
import argparse
import numpy as np
import itertools as it

import csv
import pandas as pd

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Generates RNP table with stats")
    parser.add_argument("--output_table", action="store", dest="output_table", required=True, 
                                    help="Filename of RNP output table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    #db.create_all()
    
    complexes = []
    for c in db.session.query(cdb.Complex).all():

        if c == None or not c.is_rnp():
            continue

        print c.complex_id
        complexes_dict = dict()
        complexes_dict['complex_id'] = c.complex_id
        complexes_dict['uniprot_accs'] = ' '.join([p.uniprot_acc for p in c.proteins])
        complexes_dict['genenames'] = ' '.join([p.genename() for p in c.proteins if p.genename() != None])
        complexes_dict['rna_binding_accs'] = ' '.join([p.uniprot_acc for p in c.proteins if len(p.rbp_evidences) > 0])
        complexes_dict['rna_binding_genenames'] = ' '.join([p.genename() for p in c.proteins if len(p.rbp_evidences) > 0])

        for p in c.proteins:
            print p.uniprot_acc
            print p.genename()
            print len(p.rbp_evidences) > 0
            #for rbp_e in p.rbp_evidences:
            #    print rbp_e.evidence_type

        print c.rnp_label

        complexes_dict['high_throughput'] = False
        complexes_dict['low_throughput'] = False
        complexes_dict['diffrac'] = False
        if c.rnp_label == 'ht' or c.rnp_label == 'ht_rnabinding' or c.rnp_label == 'sign_ht_rnabinding' or c.rnp_label == 'sign_ht':
            complexes_dict['high_throughput'] = True
        if c.rnp_label == 'rnabinding' or c.rnp_label == 'ht_rnabinding' or c.rnp_label == 'sign_ht_rnabinding' or c.rnp_label == 'sign_rnabinding':
            complexes_dict['low_throughput'] = True                                                                                 
        if c.rnp_label == 'sign' or c.rnp_label == 'sign20_08corr' or c.rnp_label == 'sign_ht_rnabinding' or c.rnp_label == 'sign_ht' or c.rnp_label == 'sign_rnabinding' or c.is_rnp_select():
            complexes_dict['diffrac'] = True

        complexes_dict['rnp_select'] = c.is_rnp_select()

        print ""

        complexes.append(complexes_dict)

    df = pd.DataFrame(complexes)
    df.to_csv(args.output_table, index=False)



if __name__ == "__main__":
    main()


