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


#kdrew: this is the smallest exponent before scipy.exp(x) becomes 0.0, empirically found, probably cpu specific
SMALL_EXPONENT = -745.1332
#kdrew: this is the largest exponent before scipy.exp(x) becomes inf, empirically found, probably cpu specific
LARGE_EXPONENT = 709.7827

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
    parser.add_argument("--remove_exponents", action="store_true", dest="remove_exponents", required=False, default=False,
                                    help="Do not report exponents in final output, default=False")
    parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")
    parser.add_argument("-j", "--numOfProcs", action="store", type=int, dest="numOfProcs", required=False, default=1,
                                    help="Number of processers, default=1")
    args = parser.parse_args()

    setup_log(args.logname)

    feature_table = pd.DataFrame(pd.read_csv(args.feature_matrix, sep=args.sep))
    output_df = shared_bait_feature(feature_table, args.bait_id_column, args.id_column, args.bh_correct, args.denm, args.logchoose, args.transform_pvalue, numOfProcs=args.numOfProcs, remove_exponents=args.remove_exponents)
    output_df = output_df.sort('neg_ln_pval', ascending=False)
    output_df.to_csv(args.output_file, index=False, header=True)

def shared_bait_feature(feature_table, bait_id_column, id_column, bh_correct=False, denm=False, logchoose=False, transform_pvalue=False, numOfProcs = 1, remove_exponents=False):
    print(feature_table)
    print(feature_table.columns.values)
    #kdrew: best to enforce ids as strings
    feature_table['bait_id_column_str'] = feature_table[bait_id_column].apply(str)
    print(feature_table)
    print(feature_table.columns.values)

    #kdrew: remove all rows that do not have a valid id
    q_str = "%s == %s" % (id_column, id_column)
    feature_table = feature_table.query(q_str)

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



    p = mp.Pool(numOfProcs)
    print("#### printing feature_table_geneid tables ####")
    #kdrew: generating full shared bait table for large datasets is memory intensive, break table down and do one id at a time
    gene_id_results = p.map(partial(shared_bait_feature_helper, feature_table=feature_table, id_column=id_column, ms_values=ms_values, N=N, denm=denm, logchoose=logchoose, transform_pvalue=transform_pvalue), set(feature_table[id_column].values))

    for gene_id_result in gene_id_results:
        output_dict2['gene_id1'] = output_dict2['gene_id1'] + gene_id_result['gene_id1']
        output_dict2['gene_id2'] = output_dict2['gene_id2'] + gene_id_result['gene_id2']
        output_dict2['pair_count'] = output_dict2['pair_count'] + gene_id_result['pair_count']
        output_dict2['neg_ln_pval'] = output_dict2['neg_ln_pval'] + gene_id_result['neg_ln_pval']
        output_dict2['pval'] = output_dict2['pval'] + gene_id_result['pval']
        output_dict2['exponents'] = output_dict2['exponents'] + gene_id_result['exponents']


    #kdrew: find value to transform by (ie. median exponent of all exponents) and transform_exponents for every entry
    if transform_pvalue:
        #kdrew: find min of exponents for each individual entry and then find the min of those
        min_value = np.min([np.min(ex) for ex in output_dict2['exponents']])
        print("min_value: %s" % min_value)
        #kdrew: normalize the min exponent to the SMALL_EXPONENT (-745.1332) which is the smallest exponent before scipy.exp(x) becomes 0.0
        #kdrew: this has the effect of shifting all of the small pvalues (~0.0) that evaluate to 0.0 to be within the machine precision
        transform_value = min_value - SMALL_EXPONENT 
        print("transform_value: %s" % transform_value)
                             
        #kdrew: calculate pvalue from exponents, transforming them if specified
        output_dict2['pval'] = [transform_exponents(ex, transform_value, calc_true_pval=False) for ex in output_dict2['exponents']]
        if remove_exponents:
            #kdrew: remove exponents, no longer needed
            output_dict2.pop('exponents')

        #output_dict2['neg_ln_pval'] = [-1.0*math.log(pvalue) for pvalue in output_dict2['pval']]
        output_dict2['neg_ln_pval'] = []
        for pvalue in output_dict2['pval']:
            try:
                neglpval = -1.0*math.log(pvalue)
                output_dict2['neg_ln_pval'].append(neglpval)
            except ValueError as ve:
                print(str(ve))
                output_dict2['neg_ln_pval'].append(np.nan)
            
 

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


def shared_bait_feature_helper(geneid, feature_table, id_column, ms_values, N, denm, logchoose, transform_pvalue):
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
    #ks_geneid = feature_shared_bait_table_geneid.groupby('IDs_tup')['bait_id_column_str'].nunique()
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
    return_output_dict2['exponents'] = []

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
                return_output_dict2['gene_id1'].append(str(geneid))
                return_output_dict2['gene_id2'].append(str(gene_id2))
                return_output_dict2['pair_count'].append(k)
                return_output_dict2['neg_ln_pval'].append(np.nan)
                return_output_dict2['pval'].append(np.nan)
                return_output_dict2['exponents'].append(exponents)
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

                return_output_dict2['gene_id1'].append(str(geneid))
                return_output_dict2['gene_id2'].append(str(gene_id2))
                return_output_dict2['pair_count'].append(k)
                return_output_dict2['neg_ln_pval'].append(neg_ln_p)
                return_output_dict2['pval'].append(p)
                return_output_dict2['exponents'].append(np.nan)

            except Exception as e:
                #logger.info(str(e))
                print("Exception (%s, %s): %s" % (geneid, gene_id2, str(e)))
                continue

    return return_output_dict2



if __name__ == "__main__":
    main()



