

import argparse
import itertools as it
import pandas as pd
import multiprocessing as mp
from multiprocessing import Manager

def main():

    parser = argparse.ArgumentParser(description="Tool to create table to pairwise interactions (and scores) with corresponding cluster_ids" )
    parser.add_argument("--pairwise_filename", action="store", dest="pairwise_filename", required=True,
                                            help="Filename of pairwise file, (ie. one line per ppi, third column contains score)")
    parser.add_argument("--cluster_filename", action="store", dest="cluster_filename", required=True,
                                            help="Filename of cluster file, (ie. one line per cluster)")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=False, default=None, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=',',
                                    help="Column separator for input feature matrix file, default=,")
    parser.add_argument("--header_names", action="store", nargs='+', dest="header_names", required=False, default=None,
                                    help="Names for header relating to features")
    parser.add_argument("--features", action="store", nargs='+', dest="features", required=False, default=None,
                                    help="Names of fields for features in input_feature_matrix, will test !=0.0 unless --write_value set, can also include index range (2:222)")
    parser.add_argument("--write_value", action="store_true", dest="write_value", required=False, default=False,
                                    help="Output values for features in input_feature_matrix specified by --features")
    parser.add_argument("--id_columns", action="store", nargs='+', dest="id_columns", required=False, default=None,
                                    help="Names of columns for ids in input_feature_matrix")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True,
                                    help="Output filename ")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=1,
                                    help="Number processors to use (int), default=1)")
    args = parser.parse_args()

    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())

    pairwise_interactions = dict()
    f2 = open(args.pairwise_filename, "rb")
    for line in f2.readlines():
        try:
            pairwise_interactions[frozenset([line.split()[0],line.split()[1]])] = line.split()[2]
        except IndexError:
            continue

    print "read input files"
    fout = open(args.out_filename,"wb")

    feature_table = None
    if args.feature_matrix != None:
        #kdrew: read in feature matrix and set the id columns to be strings
        feature_table = pd.read_csv(args.feature_matrix,sep=args.sep, dtype = {args.id_columns[0]:str,args.id_columns[1]:str})
        feature_table = feature_table.fillna(0.0)
        if 'frozenset_ids_str_order' not in feature_table.columns:
            #kdrew: create frozenset_ids_str_order column
            feature_table['frozenset_ids'] = map(frozenset,feature_table[[args.id_columns[0],args.id_columns[1]]].values)
            feature_table['frozenset_ids_str_order'] = feature_table['frozenset_ids'].apply(list).apply(sorted).apply(str)
            feature_table = feature_table.set_index(['frozenset_ids_str_order'])

        if args.header_names != None:
            fout.write("id1\tscore\t%s\n" % "\t".join(args.header_names))
        else:
            fout.write("id1\tscore\t%s\n" % "\t".join(args.features))
    else:
        fout.write("id1\tscore\n")

    print "Wrote Header"

    #mgr = Manager()
    #ns = mgr.Namespace()
    #ns.pairwise_interactions = pairwise_interactions
    #ns.feature_table = feature_table

    p = mp.Pool(args.procs, initializer=helper_init, initargs=({'pairwise_interactions':pairwise_interactions, 'feature_table':feature_table},))
    input_list = []
    for clustid, cluster in enumerate(clusters):
        #feature_table_trim = feature_table.query("id1 in @cluster or id2 in @cluster")
        #print feature_table_trim
        #input_list.append({'clustid':clustid, 'cluster':cluster, 'pairwise_interactions':pairwise_interactions, 'feature_table':feature_table_trim, 'features':args.features, 'write_value':args.write_value})
        input_list.append({'clustid':clustid, 'cluster':cluster, 'features':args.features, 'write_value':args.write_value})

    for cluster_out in p.imap(pairwise_interaction_helper, input_list):
        fout.write(cluster_out)

    fout.close()

def helper_init(initargs):
    #print initargs
    #kdrew: still not working properly
    init_dict = initargs
    #kdrew: global variables? yuck!
    global helper_pairwise_interactions
    global helper_feature_table
    helper_pairwise_interactions = init_dict['pairwise_interactions']
    helper_feature_table = init_dict['feature_table']


def pairwise_interaction_helper(parameter_dict):

    global helper_pairwise_interactions
    global helper_feature_table

    clustid = parameter_dict['clustid']
    cluster = parameter_dict['cluster']
    features = parameter_dict['features']
    write_value = parameter_dict['write_value']

    #pairwise_interactions = parameter_dict['pairwise_interactions']
    #feature_table = parameter_dict['feature_table']
    feature_table = helper_feature_table.query("id1 in @cluster or id2 in @cluster") 

    out_result = ""
    for prot_pair in it.combinations(cluster,2):
        id1 = "%s_%s" % (clustid, prot_pair[0])
        id2 = "%s_%s" % (clustid, prot_pair[1])
        try:
            score = helper_pairwise_interactions[frozenset([prot_pair[0],prot_pair[1]])]
        except KeyError:
            continue
        if feature_table is not None:
            out_list = []
            ids_str_order = "%s" % sorted([prot_pair[0],prot_pair[1]])
            for field in features:
                if ':' in field:
                    r1 = int(field.split(':')[0])
                    r2 = int(field.split(':')[1])
                    out_list.append(str( (feature_table.loc[ids_str_order][feature_table.columns.values[r1:r2]] != 0.0).any() ) )
                elif write_value:
                    #out_list.append(str(feature_table[feature_table['frozenset_ids_str_order'] == ids_str_order][field].values[0]))
                    out_list.append(str( feature_table.loc[ids_str_order][field] ) )
                else:
                    #kdrew: seems inefficent 
                    #out_list.append(str(feature_table[feature_table['frozenset_ids_str_order'] == ids_str_order][field].values[0] != 0.0))
                    out_list.append(str( feature_table.loc[ids_str_order][field] != 0.0))

            #fout.write("%s (pp) %s\t%s\t%s" % (id1, id2, score, "\t".join(out_list)))
            out_result = out_result + "%s (pp) %s\t%s\t%s\n" % (id1, id2, score, "\t".join(out_list))
        else:
            #fout.write("%s (pp) %s\t%s" % (id1, id2, score))
            out_result = out_result + "%s (pp) %s\t%s\n" % (id1, id2, score)

    return out_result

if __name__ == "__main__":
    main()
