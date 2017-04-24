#import sys
import argparse
import itertools as it
import pandas as pd
#sys.path.append('/project/cmcwhite/protein_complex_maps/protein_complex_maps')
#import protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Tool to create table protein ids with cluster ids, uniprot acc and gene names" )
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Output filename ")
    parser.add_argument("--orthology_filename", action="store", type=str, dest="orthology_filename", required=True,
                                            help="EggNOG orthology output file")

    parser.add_argument("--egg_filename", action="store", type=str, dest="egg_filename", required=True,
                                            help="EggNOG orthology output file")

    parser.add_argument("--arath_filename", action="store", type=str, dest="arath_filename", required=True,
                                            help="arabidopsis proteome annotations")


    #parser.add_argument("--from_id", action="store", dest="from_id", required=True, 
                                            #help="input id type, (P_ENTREZGENEID, ENSEMBL_ID,etc)")
    #parser.add_argument("--reviewed", action="store_true", dest="reviewed", required=False, default=False,
                                            #help="map only to reviewed ids, default=False")
    args = parser.parse_args()

    protid_set = set()

    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())
        map( protid_set.add, line.split() )

    #print "cluster2node_table: convert to ACC"
    #inputID2ACC_map = pu.map_protein_ids(list(protid_set), args.from_id, "ACC", reviewed=args.reviewed)
    #print "flatten_list"
    #flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
    #print "cluster2node_table: get genenames"
    #genename_map = pu.get_from_uniprot(flatten_list, 'genes')
    #print "cluster2node_table: get protein names"
    #proteinname_map = pu.get_from_uniprot(flatten_list, 'protein+names')

    d = dict()
    d['clustid'] = []
    d['GroupID'] = []
    d['clustid_key'] = []
    #d['acc'] = []
    #d['genename'] = []
    #d['proteinname'] = []
    #d['uniprot_link'] = []

    for clustid, cluster in enumerate(clusters):
        for group_id in cluster:
            clust_id_key = "%s_%s" % (clustid, group_id)

            d['clustid'].append(clustid)
            d['GroupID'].append(str(group_id))
            d['clustid_key'].append(clust_id_key)
#            try:
#                d['acc'].append(inputID2ACC_map[prot_id][0])
#                d['uniprot_link'].append("http://www.uniprot.org/uniprot/%s" % inputID2ACC_map[prot_id][0])
#                d['genename'].append(genename_map[inputID2ACC_map[prot_id][0]])
#                try:
#                    d['proteinname'].append(proteinname_map[inputID2ACC_map[prot_id][0]].strip())
#                except KeyError:
#                    d['proteinname'].append(None)
#            except IndexError:
#                d['genename'].append(None)
#                d['proteinname'].append(None)
#                d['acc'].append(None)
#                d['uniprot_link'].append(None)

    df = pd.DataFrame(d)
    ortho=pd.read_csv(args.orthology_filename, index_col=False, sep="\t")

    ortho= ortho[['GroupID', 'ProteinID', 'Annotation']]
    all_prots = ortho.groupby(['GroupID',  'Annotation'])['ProteinID'].apply(lambda x: ' '.join(x)).reset_index()
    all_prots.columns= ['GroupID', 'Eggnog_annotation', 'AllMembers']


    #print all_prots
    one_prot =  ortho.groupby(['GroupID',  'Annotation']).head(1)
    one_prot = one_prot[['GroupID', 'ProteinID']]



    all_prots = all_prots.set_index(["GroupID"])


    one_prot = one_prot.set_index(["GroupID"])

    #print all_prots
    #print one_prot


    full_table = all_prots.join(one_prot, how = "left")
    split_table = full_table['ProteinID'].str.split("|", return_type='frame')
    split_table=split_table[[1]]


    
    final_annot = full_table.join(split_table, how="left")
    print final_annot
    final_annot = final_annot[['AllMembers', 1]]
    final_annot.columns= ['AllMembers', 'Entry']
#    final_annot = final_annot.set_index(["GroupID"])
    df = df.set_index(["GroupID"])
  #  print final_annot
  #  print df

    annot1 = df.join(final_annot, how="left")
    print annot1
    egg=pd.read_csv(args.egg_filename, index_col=False, sep=",")
 

    egg = egg.set_index(["GroupID"])

    annot2 = annot1.join(egg)
    annot2 = annot2.reset_index()

    annot2 = annot2.set_index(['Entry'])

    arath=pd.read_csv(args.arath_filename, index_col=False, sep="\t")

    arath = arath.set_index(['Entry'])

    annot3 = annot2.join(arath, how="left")

    print annot3
    #annot2.reset_index()
#    egg= egg[['GroupID', 'AllMembers']]
#
#    egg = egg.set_index(['GroupID'])
#
#    prof = pd.read_csv(profile, index_col=False, sep=",")
#    prof = prof.set_index(['GroupID'])
#
#    annotated = egg.join(prof, how="right")
#    print annotated

    annot3.to_csv(args.output_filename)


if __name__ == "__main__":
    main()
