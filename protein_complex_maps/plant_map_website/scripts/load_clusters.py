from __future__ import print_function
import argparse
import plant_complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--cluster_table", action="store", dest="cluster_table", required=True, 
                                    help="Filename cluster table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()


    with(open(args.cluster_table,"rb")) as cluster_table:

        count = 1


        with open("static/data/clustid_key.csv", "w") as outfile:
            outfile.write("OrthogroupID,clustid,clustid_set,order\n")

            for line in cluster_table.readlines():
    
                if count % 100 == 0:
                    print(count)
                count = count + 1
        
                # CDM header
                # ID, cut_0.99, cut_0.95, cut_0.9, cut_0.85, cut_0.8, cut_0.75, cut_0.7, cut_0.65, cut_0.6, cut_0.55, cut_0.5, cut_0.45, cut_0.4, cut_0.35, cut_0.3, cut_0.25, cut_0.2, cut_0.15, cut_0.1
                # Select cut_0.7  cut_0.6  cut_0.45 cut_0.25cut_0.7  cut_0.6  cut_0.45 cut_0.25
                if count == 1: # Skip the header line
                    continue

                split_line = line.split(',')
    
                OrthogroupID = split_line[0]
                clustid_1 = split_line[7]
                clustid_2 = split_line[9]
                clustid_3 = split_line[12]
                clustid_4 = split_line[16]
                order = split_line[-1].strip()
               
         
                h1 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_1, clustid_set = 'clustid_1') 
                h2 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_2, clustid_set = 'clustid_2') 
                h3 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_3, clustid_set = 'clustid_3') 
                h4 = cdb.get_or_create(db, cdb.Hiercomplex, clustid = clustid_4, clustid_set = 'clustid_4') 
        
                o = cdb.get_or_create(db, cdb.Orthogroup, OrthogroupID = OrthogroupID)
                db.session.add(h1)
                db.session.add(h2)
                db.session.add(h3)
                db.session.add(h4)
        
                db.session.add(o)
                db.session.commit()
        
                ocm1 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h1.id)
                ocm2 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h2.id)
                ocm3 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h3.id)
                ocm4 = cdb.get_or_create(db, cdb.OrthogroupComplexMapping, orthogroup_key=o.id, hiercomplex_key=h4.id)
         
        
                db.session.add(ocm1)
                db.session.add(ocm2)
                db.session.add(ocm3)
                db.session.add(ocm4)    
                db.session.commit()
        
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_1, "clustid_1", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_2, "clustid_2", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_3, "clustid_3", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_4, "clustid_4", order))


if __name__ == "__main__":
    main()


