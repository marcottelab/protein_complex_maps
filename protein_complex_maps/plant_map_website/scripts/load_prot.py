from __future__ import print_function
import argparse
import time
import plant_complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads a conversion table")
    parser.add_argument("--conversion_file", action="store", dest="conversion_file", required=True,                                                      help="At least 3 column comma-separated file of eggnogID Spec ProteinID")            

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()
    count = 1
    t0  = time.time()
    with(open(args.conversion_file,"rb")) as conversion_table:
       
        ps = []
        for line in conversion_table.readlines():
            count = count + 1 
            if "ProteinID" in line:
                continue
            split_line = line.split(',')   
            OrthogroupID = split_line[0]
            ProteinID = split_line[1]
            Spec = split_line[2].strip('\n')
            o = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID).first()
            #o = cdb.get_or_create(db, cdb.Protein, ProteinID = ProteinID)
            p = cdb.Protein(ProteinID = ProteinID, Spec = Spec, OrthogroupID_key = o.id) # Plain load
            ps.append(p)  
            if count % 10000 == 0:
                 
                print(count, str(time.time() - t0))
                t0 = time.time()
                db.session.add_all(ps)
                db.session.flush()  # Gets each object updated with its .id
                db.session.commit()
                ps = []    
 
    #db.session.add_all(ps)
    #db.session.flush()  # Gets each object updated with its .id
    print("added all ProteinIDs!")
    #db.session.commit()
  

if __name__ == "__main__":

    main()
