from __future__ import print_function
import argparse
import time
import plant_complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads a score table")

    parser.add_argument("--score_file", action="store", dest="score_file", required=True,
                                    help="3 column tab-separated file of eggnogID eggnogID score")


    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    score_table = open(args.score_file,"rb")
    count = 1
    t0 = time.time()
    scrs = []
    with(open(args.score_file,"rb")) as score_table:

        for line in score_table.readlines():
            if "P_1" in line: #header
                 continue
            count = count + 1

            split_line = line.split(',')
            OrthogroupID1 = split_line[0]
            OrthogroupID2 = split_line[1]
            ScoreVal = split_line[2].strip("\n")


            o1 = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID1).first()
            o2 = db.session.query(cdb.Orthogroup).filter_by(OrthogroupID = OrthogroupID2).first()

            scr1 = cdb.get_or_create(db, cdb.Score, InteractionID = count, OrthogroupID_key=o1.id, ScoreVal = ScoreVal)
            scrs.append(scr1)

            scr2 = cdb.get_or_create(db, cdb.Score, InteractionID = count, OrthogroupID_key=o2.id, ScoreVal = ScoreVal)
            scrs.append(scr2)

            if count % 1000 == 0:
                print(count, str(time.time() - t0))
                t0 = time.time()
                db.session.add_all(scrs)
                db.session.flush()
                db.session.commit()
                scrs = []
      
    db.session.add_all(scrs)
    db.session.flush()
    db.session.commit()
    print("added all scores")

if __name__ == "__main__":
    main()
