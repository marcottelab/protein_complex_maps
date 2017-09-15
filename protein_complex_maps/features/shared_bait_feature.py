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

#kdrew: calculate pvalue using hypergeometric distribution
#kdrew: there are a few different implementations to calculate this
#kdrew: adhoc is Andrew Dalke's method of choose
#kdrew: denm does not calculate the normalization choose(N,m) but rather normalizes based on m (not a true pval)
#kdrew: logchoose calculates choose using the log transform (specifically the loggamma function)
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
            except ValueError as ve:
                print(str(ve))
                return np.nan
            r3 = logchoose_func(N,m)

            pi = scipy.exp(r1 + r2 - r3)
        else:
            pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) / misc.comb(N,m)
        pv += pi
    return pv

#kdrew: for large values the pvals will all be 0 because of machine precision and we won't have discriminating power between entries
#kdrew: returning the exponents and staying in logspace gives us a way to compare and rank entries
def hypergeometric_exponents(k,n,m,N): 
    exponents = []
    for i in range(k,int(min(m,n)+1)):
        r1 = logchoose_func(n, i)
        try:
            r2 = logchoose_func(N-n, m-i)
        except ValueError as ve:
            print(str(ve))
            r2 = np.nan
        r3 = logchoose_func(N,m)

        pi = r1 + r2 - r3
        exponents.append(pi)
    return exponents

#kdrew: transforms pval calculation to stay within machine precision
#kdrew: transform takes the form of 1/e^transform_val  * e^x (ie. implemented as x - transform_val)
#kdrew: calc_true_pval will retransform back into a true pvalue, otherwise keep in transformed space
def transform_pval(k,n,m,N, transform_val, calc_true_pval=True):
    exponents = hypergeometric_exponents(k,n,m,N)
    return transform_exponents(exponents, transform_val, calc_true_pval)

def transform_exponents(exponents, transform_val, calc_true_pval=True):
    transformed_exponents = [ex-transform_val for ex in exponents]
    transformed_pval = sum([scipy.exp(ex) for ex in transformed_exponents])
    if calc_true_pval:
        transformed_pval = transformed_pval * scipy.exp(transform_val)
    return transformed_pval


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
    parser.add_argument("--transform_pvalue", action="store_true", dest="transform_pvalue", required=False, default=False,
                                    help="Report a transformed pvalue, default=False")
    parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")
    args = parser.parse_args()

    setup_log(args.logname)

    feature_table = pd.DataFrame(pd.read_csv(args.feature_matrix, sep=args.sep))
    output_df = shared_bait_feature(feature_table, args.bait_id_column, args.id_column, args.bh_correct, args.denm, args.logchoose, args.transform_pvalue)
    output_df = output_df.sort('neg_ln_pval', ascending=False)
    output_df.to_csv(args.output_file, index=False, header=True)

def shared_bait_feature(feature_table, bait_id_column, id_column, bh_correct=False, denm=False, logchoose=False, transform_pvalue=False):
    print(feature_table)
    print(feature_table.columns.values)
    #kdrew: best to enforce ids as strings
    feature_table['bait_id_column_str'] = feature_table[bait_id_column].apply(str)
    print(feature_table)
    print(feature_table.columns.values)

    N = feature_table['bait_id_column_str'].nunique()
    feature_table['gene_id_str'] = feature_table[id_column].apply(str)
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
    output_dict2['exponents'] = []

    print("#### printing feature_table_geneid tables ####")
    #kdrew: generating full shared bait table for large datasets is memory intensive, break table down and do one id at a time
    for geneid in set(feature_table[id_column].values):
        print("geneid: %s" % geneid)
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

        #kdrew: calculate the number of experiments that each pair of proteins share
        feature_shared_bait_table_geneid = feature_shared_bait_table_geneid.reset_index()
        #ks_geneid = feature_shared_bait_table_geneid.groupby('IDs_tup')['bait_id_column_str'].nunique()
        ks_geneid = feature_shared_bait_table_geneid.groupby('gene_id2_str')['bait_id_column_str'].nunique()
        print("ks_geneid")
        print(ks_geneid)

        for gene_id2 in ks_geneid.index:
            k = ks_geneid[str(gene_id2)]
            m = ms_values[str(gene_id2)]
            n = ms_values[str(geneid)]

            #kdrew: do not calculate for the same gene
            if str(geneid) == str(gene_id2):
                continue

            if transform_pvalue:
                try:
                    exponents = hypergeometric_exponents(k,n,m,N)
                    output_dict2['gene_id1'].append(str(geneid))
                    output_dict2['gene_id2'].append(str(gene_id2))
                    output_dict2['pair_count'].append(k)
                    output_dict2['neg_ln_pval'].append(np.nan)
                    output_dict2['pval'].append(np.nan)
                    output_dict2['exponents'].append(exponents)
                except Exception as e:
                    print("Exception (%s, %s): %s" % (geneid, gene_id2, str(e)))
                    continue

            else:
                try:
                    p = pval(k,n,m,N, denm=denm, logchoose=logchoose)
                    try:
                        neg_ln_p = -1.0*math.log(p)
                    except ValueError as ve:
                        print(str(ve))
                        neg_ln_p = np.nan

                    #print("calculated pval")
                    #print("%s,%s k:%s m:%s n:%s N:%s -ln(p):%s" % (geneid, gene_id2, k, m, n, N, neg_ln_p))

                    output_dict2['gene_id1'].append(str(geneid))
                    output_dict2['gene_id2'].append(str(gene_id2))
                    output_dict2['pair_count'].append(k)
                    output_dict2['neg_ln_pval'].append(neg_ln_p)
                    output_dict2['pval'].append(p)
                    output_dict2['exponents'].append(np.nan)

                except Exception as e:
                    #logger.info(str(e))
                    print("Exception (%s, %s): %s" % (geneid, gene_id2, str(e)))
                    continue


    #kdrew: find value to transform by (ie. max exponent of all exponents) and transform_exponents for every entry
    if transform_pvalue:
        #kdrew: find max of exponents for each individual entry and then find the max of those
        max_value = max([max(ex) for ex in output_dict2['exponents']])
        print("max_value")
        print(max_value)
                             
        output_dict2['pval'] = [transform_exponents(ex, max_value, calc_true_pval=False) for ex in output_dict2['exponents']]
        output_dict2['neg_ln_pval'] = [-1.0*math.log(pval) for pval in output_dict2['pval']]
 

    if bh_correct:
        #kdrew: use Benjamini Hochberg from R to correct for multiple hypothesis testing
        stats = importr('stats')
        p_adjust = stats.p_adjust(FloatVector(output_dict2['pval']), method = 'BH')
        output_dict2['pval_corr'] = p_adjust
        output_dict2['neg_ln_pval_corr'] = []
        #kdrew: original code without error handling [-1.0*math.log(p) for p in p_adjust]
        for p in p_adjust:
            try:
                output_dict2['neg_ln_pval_corr'].append(-1.0*math.log(p))
            except ValueError as ve:
                print(str(ve))
                output_dict2['neg_ln_pval_corr'].append(np.nan)

    output_df = pd.DataFrame(output_dict2)
    return output_df

if __name__ == "__main__":
    main()



