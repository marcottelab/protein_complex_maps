import pandas as pd
import argparse


def get_groups(input_list, conversion_tbl):

    groups = []
    notfound = []

    
    for ID in input_list:
        output = conversion_tbl[conversion_tbl['ProteinID'].str.contains(ID)]
        groupID = output['ID'].tolist()[0]
        groups.append(groupID)
    
    return groups    




if __name__ == '__main__':

    input_list=["WRK58_ARATH", "Q9SR92"]
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    get_groups(input_list, conversion_tbl)
