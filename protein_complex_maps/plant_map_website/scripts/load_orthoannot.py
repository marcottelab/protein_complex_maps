import argparse
import plant_complex_db as cdb
import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Loads a annotation table")
    parser.add_argument("--annotation_file", action="store", dest="annotation_file", required=True,
                                    help="At least 3 column tab-separated file of eggnogID Arath_genename, orthoannotation")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    annotation_table = open(args.annotation_file,"rb")
    complex_id = None
    count = 1

#ID      arath_genenames arath_Entries   arath_Entry_names       arath_Protein_names     disruptions     tair_disruptionslloyd2012_LOFs  arath_functions arath_misc      pathway unipathway      BioCyc  Reactome        BRENDA  kegg_pws       ec       arath_masses    arath_protein_names     arath_GO        devstages       tissues tair    araport Annotation     orysj_genenames  orysj_Entries   orysj_Entry_names       orysj_Protein_names     orysj_disruptions       orysj_functions
#orysj_misc

    with(open(args.annotation_file,"rb")) as annotation_table:
        
        for line in annotation_table.readlines():
            if "arath_genename" in line: # Header
                continue 


            if count % 100 == 0:
                print count
            count = count + 1    
            split_line = line.split('\t')
            try:
                OrthogroupID = split_line[0]
                ArathGenenames = split_line[1]
                EggnogAnnot = split_line[24]
                Tair = split_line[23].strip("\n")
                o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
                a = cdb.get_or_create(db, cdb.Orthoannot, OrthogroupID_key = o.id, EggnogAnnot = EggnogAnnot, Tair = Tair, ArathGenenames = ArathGenenames)
                db.session.add(a)
                db.session.commit()
            except Exception as e:
                 continue

if __name__ == "__main__":

    main()
