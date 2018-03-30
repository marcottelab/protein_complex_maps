

import argparse
import itertools as it
import pandas as pd

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
                                    help="Names of fields for features in input_feature_matrix, will test !=0.0 unless --write_value set")
    parser.add_argument("--write_value", action="store_true", dest="write_value", required=False, default=False,
                                    help="Output values for features in input_feature_matrix specified by --features")
    parser.add_argument("--id_columns", action="store", nargs='+', dest="id_columns", required=False, default=None,
                                    help="Names of columns for ids in input_feature_matrix")
    parser.add_argument("--input_test_ppis", action="store", dest="test_ppis_file", required=False, default=None, 
                                    help="Filename of test ppis")
    parser.add_argument("--input_train_ppis", action="store", dest="train_ppis_file", required=False, default=None, 
                                    help="Filename of train ppis")
    parser.add_argument("--output_filename", action="store", dest="out_filename", required=True,
                                            help="Output filename ")
    args = parser.parse_args()



    clusters = []
    f = open(args.cluster_filename, "rb")
    for line in f.readlines():
        clusters.append(line.split())

    pairwise_df = pd.read_csv(args.pairwise_filename, sep='\t', header=None, names=['id1','id2','score'], dtype={'id1':str,'id2':str})
    pairwise_df['frozenset_ids'] = map(frozenset,pairwise_df[['id1','id2']].values)
    pairwise_df['frozenset_ids_str_order'] = pairwise_df['frozenset_ids'].apply(list).apply(sorted).apply(str)

    cluster_entries = []
    for clustid, cluster in enumerate(clusters):
        for prot_pair in it.combinations(cluster,2):
            id1 = "%s_%s" % (clustid, prot_pair[0])
            id2 = "%s_%s" % (clustid, prot_pair[1])
            cluster_entries.append([id1,id2,prot_pair[0], prot_pair[1]])

    cluster_df = pd.DataFrame(cluster_entries, columns=['clustid1','clustid2','id1','id2'])
    cluster_df['edge_id'] = ["%s (pp) %s" % (x[0],x[1] for x in cluster_df[['clustid1','clustid2']].values]
    cluster_df['frozenset_ids'] = map(frozenset,cluster_df[['id1','id2']].values)
    cluster_df['frozenset_ids_str_order'] = cluster_df['frozenset_ids'].apply(list).apply(sorted).apply(str)
    pairwise_df = pairwise_df.merge(cluster_df, on='frozenset_ids_str_order', how='outer', suffixes=['','_cluster'])

    columns = ['edge_id','id1','id2','score']

    if args.test_ppis_file != None:
        test_df = pd.read_csv(args.test_ppis_file, sep='\t', header=None, names=['id1','id2'], dtype={'id1':str, 'id2':str})
        test_df['frozenset_ids'] = map(frozenset,test_df[['id1','id2']].values)
        test_df['frozenset_ids_str_order'] = test_df['frozenset_ids'].apply(list).apply(sorted).apply(str)
        pairwise_df = pairwise_df.merge(test_df, on='frozenset_ids_str_order', how='outer', suffixes=['','_test'])
        pairwise_df['is_test'] = pairwise_df['id1_test'] == pairwise_df['id1_test']
        columns.append('is_test')

    if args.train_ppis_file != None:
        train_df = pd.read_csv(args.train_ppis_file, sep='\t', header=None, names=['id1','id2'], dtype={'id1':str, 'id2':str})
        train_df['frozenset_ids'] = map(frozenset,train_df[['id1','id2']].values)
        train_df['frozenset_ids_str_order'] = train_df['frozenset_ids'].apply(list).apply(sorted).apply(str)
        pairwise_df = pairwise_df.merge(train_df, on='frozenset_ids_str_order', how='outer', suffixes=['','_train'])
        pairwise_df['is_train'] = pairwise_df['id1_train'] == pairwise_df['id1_train']
        columns.append('is_train')

    if args.feature_matrix != None:
        feature_table = pd.read_csv(args.feature_matrix, sep=args.sep, dtype={args.id_columns[0]:str, args.id_columns[1]:str})
        feature_table['frozenset_ids'] = map(frozenset,feature_table[[args.id_columns[0],args.id_columns[1]]].values)
        feature_table['frozenset_ids_str_order'] = feature_table['frozenset_ids'].apply(list).apply(sorted).apply(str)
        pairwise_df = pairwise_df.merge(feature_table, on='frozenset_ids_str_order', how='outer', suffixes=['','_feature'])

        for i, field in enumerate(args.features):
            if args.write_value:
                pairwise_df[args.header_names[i]] = pairwise_df[field]
                columns.append(field)
            else:
                pairwise_df[args.header_names[i]] = pairwise_df[field] != 0.0
                columns.append(args.header_names[i])

    pairwise_df.query("edge_id == edge_id")[columns].to_csv(args.out_filename, index=False)

if __name__ == "__main__":
    main()
