
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.misc as misc
import math

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def pval(k,n,m,N):
    pv = 0.0
    for i in range(k,min(m,n)+1):
        pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) / misc.comb(N,m)
        pv += pi
    return pv


def main():
    parser = argparse.ArgumentParser(description="Calculates pval for prey prey interactions in bait prey data")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output matrix")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Separator for reading csv, default=$")
    parser.add_argument("--id_column", action="store", dest="id_column", required=False, default='gene_id',
                                    help="Name of column that specify ids in feature matrix, default=gene_id")
    parser.add_argument("--bait_id_column", action="store", dest="bait_id_column", required=False, default='bait_geneid',
                                    help="Name of column that specify bait ids in feature matrix, default=bait_geneid")

    args = parser.parse_args()

    #bioplex_feature_table = pd.read_csv("/home/kdrew/data/bioplex/150408_CDF_STAR_GRAPH_Ver2594_dollarsign.txt.test",sep='$')
    bioplex_feature_table = pd.read_csv(args.feature_matrix, sep=args.sep)


    #kdrew: merge table to itself to get pairs of proteins with same bait
    bioplex_feature_shared_bait_table = bioplex_feature_table.merge(bioplex_feature_table, on=args.bait_id_column)

    #print bioplex_feature_shared_bait_table[['gene_id_x','gene_id_y']]

    #kdrew: create set of id pairs (need set because merge generates duplicate gene pairs, also deals with order)
    df_tmp = bioplex_feature_shared_bait_table[[args.id_column+'_x',args.id_column+'_y']].apply(frozenset,axis=1)
    bioplex_feature_shared_bait_table['gene_id_set'] = df_tmp
    #bioplex_feature_shared_bait_table['gene_id_set'] = bioplex_feature_shared_bait_table[[args.id_column+'_x',args.id_column+'_y']].apply(frozenset,axis=1)[args.id_column+'_x']

    #bioplex_feature_shared_bait_table['gene_id_set'] = bioplex_feature_shared_bait_table[['gene_id_x','gene_id_y']].apply(frozenset,axis=1)['gene_id_x']
    #kdrew: generate tuple of set so groupby works, apparently cannot use cmp on sets
    bioplex_feature_shared_bait_table['gene_id_tup'] = bioplex_feature_shared_bait_table['gene_id_set'].apply(tuple)

    #kdrew: number of times pair is found (unique baits), 'k' in Hart etal 2007 
    ks = bioplex_feature_shared_bait_table.groupby('gene_id_tup').bait_geneid.nunique()
    #kdrew: number of times individual id is found (unique baits), 'm' and 'n' in Hart etal 2007 
    ms = bioplex_feature_shared_bait_table.groupby(args.id_column+'_x').bait_geneid.nunique()
    #kdrew: number of total experiments (unique baits), 'N' in Hart etal 2007 
    N = bioplex_feature_shared_bait_table.bait_geneid.nunique()


    #for gene_ids in bioplex_feature_shared_bait_table.gene_id_tup:
    output_dict = dict()
    output_dict['gene_id1'] = []
    output_dict['gene_id2'] = []
    output_dict['pair_count'] = []
    output_dict['neg_ln_pval'] = []
    for gene_ids in ks.index:
        if len(gene_ids) == 2:
            k = ks[gene_ids]
            m = ms[gene_ids[0]]
            n = ms[gene_ids[1]]
            #p = stats.hypergeom.cdf(k, N, m, n)
            p = pval(k,n,m,N)
            neg_ln_p = -1.0*math.log(p)
            #print "%s k:%s n:%s m:%s -ln(p):%s" % (gene_ids, k, m, n, neg_ln_p)

            output_dict['gene_id1'].append(gene_ids[0])
            output_dict['gene_id2'].append(gene_ids[1])
            output_dict['pair_count'].append(k)
            output_dict['neg_ln_pval'].append(neg_ln_p)


    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(args.output_file)

if __name__ == "__main__":
    main()



