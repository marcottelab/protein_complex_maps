from __future__ import print_function
import argparse

def main():

    parser = argparse.ArgumentParser(description="Loads sql tables from input files")
    parser.add_argument("--cluster_table", action="store", dest="cluster_table", required=True, 
                                    help="Filename cluster table")

    args = parser.parse_args()

    with(open(args.cluster_table,"rb")) as cluster_table:

       count = 1

       with open("clustid_key.csv", "w") as outfile:
           outfile.write("OrthogroupID,clustid,clustid_set,order\n")

    
           for line in cluster_table.readlines():

    
                if count % 100 == 0:
                    print(count)
    
                count = count + 1
    
    
                if "OrthogroupID" in line: # Skip the header line
                    continue
    
                split_line = line.split(',')
    
    
                # print line8, 11, 13, 17
                OrthogroupID = split_line[0]
                clustid_1 = split_line[7]
                clustid_2 = split_line[10]
                clustid_3 = split_line[12]
                clustid_4 = split_line[16]
                order = split_line[-1].strip()
    
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_1, "clustid_1", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_2, "clustid_2", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_3, "clustid_3", order))
                outfile.write("{},{},{},{}\n".format(OrthogroupID, clustid_4, "clustid_4", order))


if __name__ == "__main__":
    main()


