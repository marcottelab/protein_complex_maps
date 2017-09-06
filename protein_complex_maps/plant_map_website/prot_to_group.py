import pandas as pd
import argparse


def prot2group(input_list, conversion_tbl):

    df_prot = pd.DataFrame()
    for protID in input_list:
        outputprot = conversion_tbl[conversion_tbl['ProteinID'].str.contains(protID)]
        outputprot = outputprot[['ProteinID', 'ID', 'Species']]
        outputprot.columns = ['ProteinID', 'OrthogroupID', 'SpeciesCode']
        df_prot = df_prot.append(outputprot)

    #Return identifiers that are in the orthology mapping
    df_prot_html = df_prot.to_html(classes='ProtTbl', index=False)
    return df_prot_html




if __name__ == '__main__':

    input_list=["WRK58_ARATH", "Q9SR92", "W5B5G3", "BLAAA"]
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    prot2group(input_list, conversion_tbl)
