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
    return df_prot_html, df_prot['OrthogroupID']


def group2prot(input_list, conversion_tbl, spec):

    df_group = pd.DataFrame()
    for groupID in input_list:
        #print "group 2 prot groupID" + groupID
        outputgroup = conversion_tbl[conversion_tbl['ID'] == groupID]
       
        outputgroup = outputgroup[['ID', 'ProteinID', 'Species']]
        outputgroup.columns = ['OrthogroupID', 'ProteinID', 'SpeciesCode']
        df_group = df_group.append(outputgroup)

    df_group= df_group[df_group['SpeciesCode']==spec]
   
    #Return identifiers that are in the orthology mapping
    df_group_html = df_group.to_html(classes='GroupTbl', index=False)
    return df_group_html

def prot2tair(input_list, conversion_tbl):
    #Bare skeleton of a function, not written or tests
    df_prot = pd.DataFrame()
    for protID in input_list:
        outputprot = conversion_tbl[conversion_tbl['ProteinID'].str.contains(protID)]
        outputprot = outputprot[['ProteinID', 'TAIR']]
        outputprot.columns = ['ProteinID', 'TAIR']
        df_prot = df_prot.append(outputprot)

    #Return identifiers that are in the orthology mapping
    df_prot_html = df_prot.to_html(classes='ProtTbl', index=False)
    return df_prot_html, df_prot['OrthogroupID']




def make_conversion_tables(input_list, conversion_tbl, spec):
     df_prot_html, IDs = prot2group(input_list, conversion_tbl)
     df_group_html= group2prot(IDs, conversion_tbl, spec)
     return df_prot_html, df_group_html



if __name__ == '__main__':

    input_list=["WRK58_ARATH", "Q9SR92", "W5B5G3", "BLAAA"]
    conversion_tbl = pd.read_csv("all_tophits_grouplength.txt", sep="\t")

    group2group(input_list, conversion_tbl, spec)
