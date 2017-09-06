import pandas as pd
import argparse


def valid_query(input_list, conversion_tbl):

    found = []
    notfound = []

    for ID in input_list:
        output = conversion_tbl[conversion_tbl['ProteinID'].str.contains(ID)]

        if len(output) == 1:
           found.append(ID)
        elif len(output) == 0:
           notfound.append(ID)
        

    if len(notfound) > 0:
        print ",".join(notfound) +" not found"
    #Return identifiers that are in the orthology mapping

    return found 




if __name__ == '__main__':

    input_list=["WRK58_ARATH", "Q9SR92", "W5B5G3", "BLAAA"]
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    valid_query(input_list, conversion_tbl)
