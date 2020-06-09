from __future__ import print_function
import sys
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.misc as misc
import scipy
import math
import logging 
import multiprocessing as mp
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from functools import partial

import mpmath as mpm


pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def pval(k,n,m,N):
    pv = 0.0
    for i in range(k,int(min(m,n)+1)):
        pi = ( mpm.binomial(n,i) * mpm.binomial((N-n), (m-i)) ) / mpm.binomial(N,m)
        pv += pi
    return pv

#kdrew: calculate pvalue using hypergeometric distribution
#kdrew: there are a few different implementations to calculate this
#kdrew: adhoc is Andrew Dalke's method of choose
#kdrew: denm does not calculate the normalization choose(N,m) but rather normalizes based on m (not a true pval)
#kdrew: logchoose calculates choose using the log transform (specifically the loggamma function)
def pval_old(k,n,m,N, adhoc=False, denm=False, logchoose=False):
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
            except ValueError as ve:
                print(str(ve))
                return np.nan
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
    parser.add_argument("--abundance_column", action="store", dest="abundance_column", required=False, default='abundance',
                                    help="Name of column that specifies abundance in feature matrix, default=abundance")
    parser.add_argument("--abundance_threshold", action="store", type=float, dest="abundance_threshold", required=False, default=0,
                                    help="Threshold abundance above (>=) a given value, default=0")
    parser.add_argument("--bh_correct", action="store_true", dest="bh_correct", required=False, default=False,
                                    help="Benjamini-Hochberg correct pvals, default=False")
    parser.add_argument("--use_abundance", action="store_true", dest="use_abundance", required=False, default=False,
                                    help="Use abundance measures when calculating hypergeometric test, default=False (presence-absence)")
    parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")
    parser.add_argument("-j", "--numOfProcs", action="store", type=int, dest="numOfProcs", required=False, default=1,
                                    help="Number of processers, default=1")
    args = parser.parse_args()

    setup_log(args.logname)

    feature_table = pd.DataFrame(pd.read_csv(args.feature_matrix, sep=args.sep, converters={args.id_column: str, args.bait_id_column: str }))
    print("Thresholding %s >= %s" % (args.abundance_column, args.abundance_threshold))
    feature_table = feature_table[feature_table[args.abundance_column] >= args.abundance_threshold]
    output_df = shared_bait_feature(feature_table, args.bait_id_column, args.id_column, args.abundance_column, args.bh_correct, args.use_abundance, numOfProcs=args.numOfProcs)
    output_df = output_df.sort('neg_ln_pval', ascending=False)
    output_df.to_csv(args.output_file, index=False, header=True)

def shared_bait_feature(feature_table, bait_id_column, id_column, abundance_column='abundance', bh_correct=False, use_abundance=False, numOfProcs = 1):
    print(feature_table)
    print(feature_table.columns.values)
    #kdrew: best to enforce ids as strings
    feature_table['bait_id_column_str'] = feature_table[bait_id_column].apply(str)
    print(feature_table)
    print(feature_table.columns.values)


    #kdrew: remove all rows that do not have a valid id
    q_str = "%s == %s" % (id_column, id_column)
    feature_table = feature_table.query(q_str)

    feature_table['gene_id_str'] = feature_table[id_column].apply(str)

    if use_abundance:
        #kdrew: convert all abundances to ints
        feature_table['abundance_int'] = feature_table[abundance_column].apply(int)

        #kdrew: still working out details about what N should be, whether max abundance of every experiment or sum abundance of every experiment
        #N = feature_table.groupby("bait_id_column_str")['abundance_int'].sum().sum()
        N = feature_table.groupby("bait_id_column_str")['abundance_int'].max().sum()

        #kdrew: number of experiments each gene_id is present in
        ms_values = feature_table.groupby('gene_id_str')['abundance_int'].sum()

    else:
        N = feature_table['bait_id_column_str'].nunique()

        #kdrew: number of experiments each gene_id is present in
        ms_values = feature_table.groupby('gene_id_str')[bait_id_column].nunique()

    print("ms_values")
    print(ms_values)
    #feature_table = feature_table.drop(bait_id_column)
    feature_table = feature_table.set_index('bait_id_column_str')

    #kdrew: create dictionary to store output, gets converted into pandas dataframe at the end
    output_dict2 = dict()
    output_dict2['gene_id1'] = []
    output_dict2['gene_id2'] = []
    output_dict2['pair_count'] = []
    output_dict2['neg_ln_pval'] = []
    output_dict2['pval'] = []

    p = mp.Pool(numOfProcs)
    print("#### printing feature_table_geneid tables ####")
    #kdrew: generating full shared bait table for large datasets is memory intensive, break table down and do one id at a time
    gene_id_results = p.map(partial(shared_bait_feature_helper, feature_table=feature_table, id_column=id_column, use_abundance=use_abundance, ms_values=ms_values, N=N), set(feature_table[id_column].values))

    for gene_id_result in gene_id_results:
        output_dict2['gene_id1'] = output_dict2['gene_id1'] + gene_id_result['gene_id1']
        output_dict2['gene_id2'] = output_dict2['gene_id2'] + gene_id_result['gene_id2']
        output_dict2['pair_count'] = output_dict2['pair_count'] + gene_id_result['pair_count']
        output_dict2['neg_ln_pval'] = output_dict2['neg_ln_pval'] + gene_id_result['neg_ln_pval']
        output_dict2['pval'] = output_dict2['pval'] + gene_id_result['pval']

    if bh_correct:
        #kdrew: use Benjamini Hochberg from R to correct for multiple hypothesis testing
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(output_dict2['pval']), method = 'BH')
        output_dict2['pval_corr'] = p_adjust
        output_dict2['neg_ln_pval_corr'] = []
        #kdrew: original code without error handling [-1.0*math.log(p) for p in p_adjust]
        for p in p_adjust:
            try:
                output_dict2['neg_ln_pval_corr'].append(-1.0*mpm.log(p))
            except ValueError as ve:
                print(str(ve))
                output_dict2['neg_ln_pval_corr'].append(np.nan)

    output_df = pd.DataFrame(output_dict2)

    output_df['fset_ids'] = [frozenset(x) for x in output_df[['gene_id1','gene_id2']].values]
    output_df = output_df.groupby(output_df.fset_ids).first()
    output_df = output_df[['gene_id1','gene_id2','neg_ln_pval','pair_count','pval']]
    return output_df


def shared_bait_feature_helper(geneid, feature_table, id_column, use_abundance, ms_values, N ):
    print("geneid: %s" % geneid)
    sys.stdout.flush()
    #kdrew: slice feature_table to only have a single id
    feature_table_geneid = feature_table[feature_table[id_column] == geneid] 
    #kdrew: join the sliced feature table with the entire feature table on the experiment column (ie. bait_id_column_str)
    feature_shared_bait_table_geneid = feature_table_geneid.join(feature_table, rsuffix="_right")
    #print (feature_table_geneid)

    #kdrew: make ids into strings and generate tuples of pairs of genes
    feature_shared_bait_table_geneid['gene_id1_str'] = feature_shared_bait_table_geneid[id_column].apply(str)
    feature_shared_bait_table_geneid['gene_id2_str'] = feature_shared_bait_table_geneid[id_column+'_right'].apply(str)
    feature_shared_bait_table_geneid['IDs'] = map(sorted, zip(feature_shared_bait_table_geneid['gene_id1_str'].values, feature_shared_bait_table_geneid['gene_id2_str'].values))
    feature_shared_bait_table_geneid['IDs_tup'] = feature_shared_bait_table_geneid['IDs'].apply(tuple)
    print (feature_shared_bait_table_geneid)
    sys.stdout.flush()

    #kdrew: calculate the number of experiments that each pair of proteins share
    feature_shared_bait_table_geneid = feature_shared_bait_table_geneid.reset_index()

    if use_abundance:
        #kdrew: still working out proper way of calculating k, sum isn't going to work because k could end up > n which makes n choose k = 0.0
        feature_shared_bait_table_geneid['min_abundance'] = feature_shared_bait_table_geneid[['abundance_int','abundance_int_right']].min(axis=1)
        #feature_shared_bait_table_geneid['sum_abundance'] = feature_shared_bait_table_geneid[['abundance_int','abundance_int_right']].sum(axis=1)
        #ks_geneid = feature_shared_bait_table_geneid.groupby('gene_id2_str')['sum_abundance'].sum()
        ks_geneid = feature_shared_bait_table_geneid.groupby('gene_id2_str')['min_abundance'].sum()

    else:
        ks_geneid = feature_shared_bait_table_geneid.groupby('gene_id2_str')['bait_id_column_str'].nunique()

    print("ks_geneid")
    print(ks_geneid)
    sys.stdout.flush()

    return_output_dict2 = dict()
    return_output_dict2['gene_id1'] = []
    return_output_dict2['gene_id2'] = []
    return_output_dict2['pair_count'] = []
    return_output_dict2['neg_ln_pval'] = []
    return_output_dict2['pval'] = []

    for gene_id2 in ks_geneid.index:
        k = ks_geneid[str(gene_id2)]
        m = ms_values[str(gene_id2)]
        n = ms_values[str(geneid)]

        #kdrew: do not calculate for the same gene
        if str(geneid) == str(gene_id2):
            continue

        try:
            p = pval(k,n,m,N)
            try:
                neg_ln_p = -1.0*mpm.log(p)
            except ValueError as ve:
                print(str(ve))
                neg_ln_p = np.nan

            #print("calculated pval")
            #print("%s,%s k:%s m:%s n:%s N:%s -ln(p):%s" % (geneid, gene_id2, k, m, n, N, neg_ln_p))

            return_output_dict2['gene_id1'].append(str(geneid))
            return_output_dict2['gene_id2'].append(str(gene_id2))
            return_output_dict2['pair_count'].append(k)
            return_output_dict2['neg_ln_pval'].append(neg_ln_p)
            return_output_dict2['pval'].append(p)

        except Exception as e:
            #logger.info(str(e))
            print("Exception (%s, %s): %s" % (geneid, gene_id2, str(e)))
            continue

    return return_output_dict2



if __name__ == "__main__":
    main()



