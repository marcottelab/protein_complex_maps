from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.misc as misc
import scipy
import math
import logging 
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector



pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def pval(k,n,m,N, adhoc=False, denm=False, logchoose=False):
    #logger.info("%s %s %s %s" % (k,n,m,N))
    pv = 0.0
    for i in range(k,int(min(m,n)+1)):
        #kdrew: use adhoc choose method instead of scipy.misc.comb
        if adhoc:
            pi = ( choose(n,i) * choose((N-n), (m-i)) ) / choose(N,m)
        elif denm:
            pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) / m
        elif logchoose:
            r1 = logchoose_func(n, i)
            try:
                r2 = logchoose_func(N-n, m-i)
            except ValueError:
                return 0
            r3 = logchoose_func(N,m)

            pi = scipy.exp(r1 + r2 - r3)
        else:
            pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) / misc.comb(N,m)
        pv += pi
    return pv

def logchoose_func(n,k):
    """
    #kdrew: similar to this code 
    https://gist.github.com/JudoWill/1051335
    """
    if 0 <= k <= n:
        try:
            lgn1 = math.lgamma(n+1)
            lgk1 = math.lgamma(k+1)
            lgnk1 = math.lgamma(n-k+1)
        except ValueError:
            #print ni,ki
            raise ValueError

        return lgn1 - (lgnk1 + lgk1)

    else:
        return 0

#kdrew: from http://stackoverflow.com/questions/33468821/is-scipy-misc-comb-faster-than-an-ad-hoc-binomial-computation
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def setup_log(logname):

    LOG_FILENAME = logname
    global logger
    logger = logging.getLogger()
    filehandler = logging.FileHandler(LOG_FILENAME, mode="w")
    formatter = logging.Formatter('%(message)s')
    filehandler.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(filehandler)


    logger.setLevel(logging.INFO)


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
    parser.add_argument("--bh_correct", action="store_true", dest="bh_correct", required=False, default=False,
                                    help="Benjamini-Hochberg correct pvals")
    parser.add_argument("--denm", action="store_true", dest="denm", required=False, default=False,
                                    help="Replace denominator of hypergeometric calculation with m-value instead of N choose m, pvalues are nonsense with this option, default=False")
    parser.add_argument("--logchoose", action="store_true", dest="logchoose", required=False, default=False,
                                    help="Use logs to deal with large values when calculating choose, default=False")
    parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")
    args = parser.parse_args()

    setup_log(args.logname)

    feature_table = pd.DataFrame(pd.read_csv(args.feature_matrix, sep=args.sep))
    output_df = shared_bait_feature(feature_table, args.bait_id_column, args.id_column, args.bh_correct, args.denm, args.logchoose)
    output_df = output_df.sort('neg_ln_pval', ascending=False)
    output_df.to_csv(args.output_file, index=False, header=True)

def shared_bait_feature(feature_table, bait_id_column, id_column, bh_correct=False, denm=False, logchoose=False):
    print(feature_table)
    print(feature_table.columns.values)
    #kdrew: best to enforce ids as strings
    feature_table['bait_id_column_str'] = feature_table[bait_id_column].apply(str)
    print(feature_table)
    print(feature_table.columns.values)
    #feature_table = feature_table.drop(bait_id_column)
    feature_table = feature_table.set_index('bait_id_column_str')


    print("#### printing feature_table_geneid tables ####")
    #kdrew: generating full shared bait table for large datasets is memory intensive, break table down and do one id at a time
    for geneid in set(feature_table[id_column].values):
        print("geneid: %s" % geneid)
        #kdrew: slice feature_table to only have a single id
        feature_table_geneid = feature_table[feature_table[id_column] == geneid] 
        #kdrew: join the sliced feature table with the entire feature table on the experiment column (ie. bait_id_column_str)
        feature_shared_bait_table_geneid = feature_table_geneid.join(feature_table, rsuffix="_right")
        #print (feature_table_geneid)
        feature_shared_bait_table_geneid['gene_id1_str'] = feature_shared_bait_table_geneid[id_column].apply(str)
        feature_shared_bait_table_geneid['gene_id2_str'] = feature_shared_bait_table_geneid[id_column+'_right'].apply(str)
        feature_shared_bait_table_geneid['IDs'] = map(sorted, zip(feature_shared_bait_table_geneid['gene_id1_str'].values, feature_shared_bait_table_geneid['gene_id2_str'].values))
        feature_shared_bait_table_geneid['IDs_tup'] = feature_shared_bait_table_geneid['IDs'].apply(tuple)
        print (feature_shared_bait_table_geneid)
        feature_shared_bait_table_geneid = feature_shared_bait_table_geneid.reset_index()
        ks_geneid = feature_shared_bait_table_geneid.groupby('IDs_tup')['bait_id_column_str'].nunique()
        print("ks_geneid: %s" % ks_geneid)
        #print("feature_table")
        #print(feature_table)
        ms_geneid = feature_table[feature_table[id_column] == geneid].index.nunique()
        print("ms_geneid: %s" % ms_geneid)

    print("#### done printing feature_table_geneid tables ####")

    #cmcwhite: join table to itself to get pairs of proteins with same bait
    feature_shared_bait_table = feature_table.join(feature_table, rsuffix="_right")
    print(feature_shared_bait_table)
    print("$$$$ printed full feature_shared_bait_table $$$$")

    #logger.info(feature_shared_bait_table)
    #print(feature_shared_bait_table)

    feature_shared_bait_table = feature_shared_bait_table.reset_index()

    feature_shared_bait_table['gene_id1_str'] = feature_shared_bait_table[id_column].apply(str)
    feature_shared_bait_table['gene_id2_str'] = feature_shared_bait_table[id_column+'_right'].apply(str)

    feature_shared_bait_table['IDs'] = map(sorted, zip(feature_shared_bait_table['gene_id1_str'].values, feature_shared_bait_table['gene_id2_str'].values))

    print(feature_shared_bait_table)
    print(feature_shared_bait_table.columns.values)
  
    ##kdrew: generate tuple of set so groupby works, apparently cannot use cmp on sets
    feature_shared_bait_table['IDs_tup'] = feature_shared_bait_table['IDs'].apply(tuple)

    #kdrew: number of times pair is found (unique baits), 'k' in Hart etal 2007 
    ks = feature_shared_bait_table.groupby('IDs_tup')['bait_id_column_str'].nunique()
    print("KS KS KS KS KS KS KS KS")
    print(ks)
    #kdrew: number of times individual id is found (unique baits), 'm' and 'n' in Hart etal 2007 
    ms = feature_shared_bait_table.groupby('gene_id1_str')['bait_id_column_str'].nunique()
    print("MS MS MS MS MS MS MS MS")
    print(ms)
    #kdrew: number of total experiments (unique baits), 'N' in Hart etal 2007 
    N = feature_shared_bait_table['bait_id_column_str'].nunique()
    #print(ks, ms, N)
    #for gene_ids in bioplex_feature_shared_bait_table.gene_id_tup:
    output_dict = dict()
    output_dict['gene_id1'] = []
    output_dict['gene_id2'] = []
    output_dict['pair_count'] = []
    output_dict['neg_ln_pval'] = []
    output_dict['pval'] = []

    #logger.info(ks)
    #print(ks)
    for gene_ids_str in ks.index:
        #print(gene_ids_str)
        #gene_ids_clean = gene_ids_str.translate(None, "[\'],")
        #gene_ids = gene_ids_clean.split()

        gene_ids = list(gene_ids_str)
        #print(gene_ids)
        if len(gene_ids) == 2:
            #logger.info(gene_ids) #cdm print out pairs of gene IDs
            #print(gene_ids) #cdm print out pairs of gene IDs
            k = ks[gene_ids_str]
            m = ms[gene_ids[0]]
            n = ms[gene_ids[1]]
            #cdm:We don't want coelution of the same protein...
            if gene_ids[0]==gene_ids[1]:
                 continue
            #print k
            #print m
            #print n
            #p = stats.hypergeom.cdf(k, N, m, n)
            try:
               p = pval(k,n,m,N, denm=denm, logchoose=logchoose)
               neg_ln_p = -1.0*math.log(p)
               #print("%s k:%s n:%s m:%s -ln(p):%s" % (gene_ids, k, m, n, neg_ln_p))

               output_dict['gene_id1'].append(gene_ids[0])
               output_dict['gene_id2'].append(gene_ids[1])
               output_dict['pair_count'].append(k)
               output_dict['neg_ln_pval'].append(neg_ln_p)
               output_dict['pval'].append(p)
            except Exception as e:
                 #logger.info(str(e))
                 print(str(e))
                 continue
             
 

    if bh_correct:
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(output_dict['pval']), method = 'BH')
        output_dict['pval_corr'] = p_adjust
        output_dict['neg_ln_pval_corr'] = [-1.0*math.log(p) for p in p_adjust]

    #kdrew: not sure why this is here, CDM?
    #output_df= output_df[['gene_id1','gene_id2','neg_ln_pval']]
    
    output_df = pd.DataFrame(output_dict)
    return output_df

if __name__ == "__main__":
    main()



