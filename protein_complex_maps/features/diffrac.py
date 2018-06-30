
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats

from sklearn.ensemble import RandomForestClassifier
#from sklearn.mixture import GaussianMixture
from sklearn.mixture import GMM
from statsmodels.sandbox.stats.multicomp import fdrcorrection0


from pyemd import emd

import protein_complex_maps.features.ExtractFeatures.Features as eff

def main():

    parser = argparse.ArgumentParser(description="Calculate difference features between two fractionation experiments")
    parser.add_argument("--elution_files", action="store", nargs='+', dest="elution_files", required=True, 
                                    help="Elution files (.elut)")
    parser.add_argument("--features", action="store", nargs='+', dest="features", required=False, default=['diffrac'],
                                    help="Features to calculate: diffrac (L1-norm of difference) diffrac_percent diffrac_normalized pearsonr poisson mean_abundance emd zscore sliding_zscore fdr_correct sliding_fdr_correct")
    parser.add_argument("--annotated_list", action="store", dest="annotated_list", required=False, default=None, 
                                    help="Filename of annotated ids, used for calculating zscores from compliment of list, default=None")
    parser.add_argument("--contaminate_tag", action="store", dest="contaminate_tag", required=False, default='CONTAMINANT', 
                                    help="Filters entries with tag, default=CONTAMINANT")
    parser.add_argument("--use_gmm", action="store_true", dest="use_gmm", required=False, default=False, 
                                    help="Fit sliding window distributions to Gaussian Mixture Model and use largest gaussian for calculating zscore, default=False")
    parser.add_argument("--log_transform", action="store_true", dest="log_transform", required=False, default=False, 
                                    help="Use the log transform of the diffrac score to calculate sliding zscore, default=False")
    parser.add_argument("--window_size", action="store", type=int, dest="window_size", required=False, default=100, 
                                    help="Window size to use for calculating sliding zscore, default=100")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")


    args = parser.parse_args()
    
    elutions = []
    for efile in args.elution_files:
        elut = eff.Elut()
        elut.load(efile,format='tsv')
        elut.threshold(thresh=1)
        elutions.append(elut)

    feature_df = pd.DataFrame()
    if len(elutions) >= 2:

        if 'diffrac' in args.features:
            feature_series = calc_diffrac(elutions[0], elutions[1], normalize_totalCounts=False)
            feature_series.name = 'diffrac'
            feature_df = join_feature(feature_df,feature_series)
        if 'diffrac_percent' in args.features:
            feature_series = calc_diffrac(elutions[0], elutions[1], percent_totalCounts=True)
            feature_series.name = 'diffrac_percent'
            feature_df = join_feature(feature_df,feature_series)
        if 'diffrac_normalized' in args.features:
            feature_series = calc_diffrac(elutions[0], elutions[1], normalize_totalCounts=True)
            feature_series.name = 'diffrac_normalized'
            feature_df = join_feature(feature_df,feature_series)
        if 'emd' in args.features:
            feature_series = calc_emd(elutions[0], elutions[1])
            feature_series.name = 'emd'
            feature_df = join_feature(feature_df,feature_series)
        if 'pearsonr' in args.features:
            feature_series = calc_correlation(elutions[0], elutions[1], correlation_func=lambda x,y: stats.pearsonr(x,y)[0])
            feature_series.name = 'pearsonr'
            feature_df = join_feature(feature_df,feature_series)
        if 'poisson' in args.features:
            print("WARNING: poisson not implemented")
            #feature_series = calc_correlation(elutions[0], elutions[1])
            #feature_series.name = 'poisson'
            #feature_df = join_feature(feature_df,feature_series)
        if 'mean_abundance' in args.features:
            feature_series = calc_mean_abundance(elutions[0], elutions[1])
            feature_series.name = 'mean_abundance'
            feature_df = join_feature(feature_df,feature_series)

        if args.annotated_list != None:
            #kdrew: add in training labels
            annotated_df = pd.read_table(args.annotated_list, header=None, names=['annotated'])
            annotated = [i in annotated_df['annotated'].values for i in feature_df.index]
            feature_df['annotated'] = annotated

        print len(feature_df)
        feature_df = feature_df[~feature_df.index.str.contains('CONTAMINANT')]
        print len(feature_df)

        if 'zscore' in args.features:
            if 'diffrac_normalized' not in args.features:
                #kdrew: calculating diffrac_normalized
                feature_series = calc_diffrac(elutions[0], elutions[1], normalize_totalCounts=False)
                feature_series.name = 'diffrac'
                feature_df = join_feature(feature_df,feature_series)
            feature_series = calc_zscore(feature_df)
            feature_series.name = 'zscore'
            feature_df = join_feature(feature_df,feature_series)

        if 'sliding_zscore' in args.features:
            feature_series = calc_sliding_zscore(feature_df, window=args.window_size, use_gmm=args.use_gmm, log_transform=args.log_transform)
            feature_series.name = 'sliding_zscore'
            feature_df = join_feature(feature_df,feature_series)

        if 'fdr_correct' in args.features:
            fdr_df = calc_fdr_correct(feature_df)
            feature_df = join_feature(feature_df,fdr_df)

        if 'sliding_fdr_correct' in args.features:
            sliding_fdr_df = calc_sliding_fdr_correct(feature_df)
            feature_df = join_feature(feature_df, sliding_fdr_df)

        if args.out_filename != None:
            feature_df.sort_values(args.features[0], ascending=False).to_csv(args.out_filename)
        else:
            print feature_df.sort_values(args.features[0], ascending=False)

def join_feature(df,feature):
    return df.join(feature, how='outer')

def calc_diffrac(elut1, elut2, percent_totalCounts=False, normalize_totalCounts=False):

    #kdrew: set columns to be the same, do some error checking to ensure lengths match, also if any realignment is necessary this is the place to do it.
    assert(len(elut2.df.columns) == len(elut1.df.columns))
    elut2.df.columns = elut1.df.columns

    #kdrew: add empty rows for the ids in elut1 that are not in elut2 and vice versa
    elut1_ids = set(elut1.df.index)
    elut2_ids = set(elut2.df.index)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut1_not_elut2_ids = elut1_ids - elut2_ids
    elut1_not_elut2 = elut1.df.loc[list(elut1_not_elut2_ids)]
    elut1_not_elut2[:] = 0.0
    elut2.df = elut2.df.append(elut1_not_elut2)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut2_not_elut1_ids = elut2_ids - elut1_ids
    elut2_not_elut1 = elut2.df.loc[list(elut2_not_elut1_ids)]
    elut2_not_elut1[:] = 0.0
    elut1.df = elut1.df.append(elut2_not_elut1)

    elut_diff = elut1.df.subtract(elut2.df)
    diffrac_sum = np.abs(elut_diff).sum(axis='columns')

    #kdrew: measures how much of the total counts shifted, 1.0 total shift -> 0.0 no shift
    if percent_totalCounts:
        diffrac_sum = diffrac_sum/(elut1.df.sum(axis='columns') + elut2.df.sum(axis='columns'))
    elif normalize_totalCounts:
        diffrac_sum = diffrac_sum * diffrac_sum/(elut1.df.sum(axis='columns') + elut2.df.sum(axis='columns'))

    return diffrac_sum

def calc_emd(elut1, elut2): 

    #kdrew: set columns to be the same, do some error checking to ensure lengths match, also if any realignment is necessary this is the place to do it.
    assert(len(elut2.df.columns) == len(elut1.df.columns))
    elut2.df.columns = elut1.df.columns

    #kdrew: add empty rows for the ids in elut1 that are not in elut2 and vice versa
    elut1_ids = set(elut1.df.index)
    elut2_ids = set(elut2.df.index)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut1_not_elut2_ids = elut1_ids - elut2_ids
    elut1_not_elut2 = elut1.df.loc[list(elut1_not_elut2_ids)]
    elut1_not_elut2[:] = 0.0
    elut2.df = elut2.df.append(elut1_not_elut2)

    #kdrew: add rows in elut1 in elut2 as 0.0
    elut2_not_elut1_ids = elut2_ids - elut1_ids
    elut2_not_elut1 = elut2.df.loc[list(elut2_not_elut1_ids)]
    elut2_not_elut1[:] = 0.0
    elut1.df = elut1.df.append(elut2_not_elut1)

    #kdrew: setup distance matrix, every transition costs 1.0
    dmat = np.ones((len(elut2.df.columns),len(elut2.df.columns)))
    #kdrew: make identity transitions cost 0.0
    np.fill_diagonal(dmat,0.0)

    emd_results = []
    for idx in elut1.df.index:
        x = np.ascontiguousarray(elut1.df.loc[idx])
        #print(x)
        y = np.ascontiguousarray(elut2.df.loc[idx])
        #print(y)
        emd_result = emd(x, y, dmat)
        emd_results.append(emd_result)

    #kdrew: annoying trick to compare the two dataframes using a function
    #emd_result = elut1.df.apply(lambda x: emd_func(x, elut2.df, dmat), axis=1)
    emd_results = pd.Series(emd_results)
    emd_results.index = elut1.df.index
    return emd_results

def emd_func(x, df2, dmat):
    print(x.values.flags)
    y = df2.loc[x.name]
    print(y)
    emd_result = emd(x.values, y.values, dmat)
    print emd_result
    return emd_result


def calc_correlation(elut1, elut2, correlation_func=stats.pearsonr, default=0.0):
    intersection_ids = set(elut1.df.index).intersection(set(elut2.df.index))
    union_ids = set(elut1.df.index).union(set(elut2.df.index))

    correlation_dict = {uid:default for uid in union_ids}

    for uid in intersection_ids:
        pcoeff = correlation_func(elut1.df.ix[uid],elut2.df.ix[uid])
        correlation_dict[uid] = pcoeff

    df = pd.Series(correlation_dict.values(), index=correlation_dict.keys())
    return df

def calc_mean_abundance(elut1, elut2): 
    elut1_sum = elut1.df.sum(axis=1)
    elut2_sum = elut2.df.sum(axis=1)
    df = (elut1_sum.add(elut2_sum, fill_value=0.0)) / 2.0
    return df

def calc_zscore(feat_df): 
    if 'annotated' in feat_df.columns:
        mean = np.mean(feat_df.query("~annotated")['diffrac_normalized'])
        std = np.std(feat_df.query("~annotated")['diffrac_normalized'])
    else:
        print "WARNING: Couldn't find column 'annotated', using all rows for distribution"
        mean = np.mean(feat_df['diffrac_normalized'])
        std = np.std(feat_df['diffrac_normalized'])
    df = (feat_df['diffrac_normalized'] - mean)/std
    return df

#kdrew: min_weight_threshold : mixture model weight has to be above threshold in order to use
def calc_sliding_zscore(feat_df, window=100, use_gmm=False, min_weight_threshold=0.6, log_transform=False): 
    sliding_zscore_dict = dict()

    for id1 in feat_df.sort_values("mean_abundance",ascending=False).query("mean_abundance == mean_abundance").index:
        i_abnd = feat_df.ix[id1]['mean_abundance']
        #kdrew: entries greater than current id
        if 'annotated' in feat_df.columns:
            gt_entries = feat_df.query("~annotated and (mean_abundance > %s)" % i_abnd).sort_values('mean_abundance')['mean_abundance']
            lt_entries = feat_df.query("~annotated and (mean_abundance < %s)" % i_abnd).sort_values('mean_abundance', ascending=False)['mean_abundance']
        else:
            print "WARNING: Couldn't find column 'annotated', using all rows for distribution"
            gt_entries = feat_df.query("(mean_abundance >= %s)" % i_abnd).sort_values('mean_abundance')['mean_abundance']
            lt_entries = feat_df.query("(mean_abundance < %s)" % i_abnd).sort_values('mean_abundance', ascending=False)['mean_abundance']

        print "gt_entries"
        print gt_entries
        print "lt_entries"
        print lt_entries
        h = window
        j = window
        #kdrew: if not enough entries, adjust the other index
        if len(gt_entries) < h: j = j + (h - len(gt_entries)); h = len(gt_entries)
        if len(lt_entries) < j: h = h + (j - len(lt_entries)); j = len(lt_entries)

        entries = list(gt_entries.index[:h]) + list(lt_entries.index[:j])
        if log_transform:
            diffrac_normalized_list = (feat_df.ix[entries]['diffrac_normalized'].fillna(0.0)+0.1).apply(np.log10)
        else:
            diffrac_normalized_list = feat_df.ix[entries]['diffrac_normalized'].values
        if use_gmm:
            #kdrew: probably should be careful about using GMM's interface, originally was using GaussianMixture but that only exists in newer versions of sklearn
            #kdrew: create two models, one with a single gaussian and one with two gaussians
            gmm1 = GMM(n_components=1, covariance_type='spherical').fit(diffrac_normalized_list.reshape(-1,1))
            gmm2 = GMM(n_components=2, covariance_type='spherical').fit(diffrac_normalized_list.reshape(-1,1))
            #kdrew: Calculate their Baysian Information Criterion which penalizes additional parameters 
            gmm1_bic = gmm1.bic(diffrac_normalized_list.reshape(-1,1))
            gmm2_bic = gmm2.bic(diffrac_normalized_list.reshape(-1,1))
            print "gmm1 BIC: %s" % gmm1_bic
            print "gmm2 BIC: %s" % gmm2_bic

            print "gmm2.means_ %s" % gmm2.means_
            min_mean_model = np.argmin(gmm2.means_)
            print "gmm2.weights_ %s" % gmm2.weights_
            max_weight_model = np.argmax(gmm2.weights_)
            print "[np.sqrt(x) for x in gmm2.covars_] %s" % [np.sqrt(x) for x in gmm2.covars_]

            #kdrew: use Baysian Information Criterion for model selection
            #kdrew: also tests to make sure the model with the lowest mean is the dominant peak, 
            #kdrew: also checks that the dominant peak is above some threshold of dominance (might not be necessary anymore with BIC selection but nice to have option
            if gmm1_bic < gmm2_bic or min_mean_model != max_weight_model or np.max(gmm2.weights_) < min_weight_threshold:
                print "WARNING: Two-component GMM has higher BIC than one-component GMM or \n highest weighted model does not equal lowest mean model or \n min_weight_threshold not satisfied, *not* using gaussian mixture model"
                mean_tmp = np.mean(diffrac_normalized_list)
                std_tmp = np.std(diffrac_normalized_list)
            else:
                mean_tmp = gmm2.means_[min_mean_model][0]
                std_tmp = np.sqrt(gmm2.covars_[min_mean_model][0])

        else:
            mean_tmp = np.mean(diffrac_normalized_list)
            std_tmp = np.std(diffrac_normalized_list)

        print "diffrac_normalized_list %s" % diffrac_normalized_list
        print "id: %s mean_tmp: %s std_tmp: %s" % (id1, mean_tmp, std_tmp)
        if log_transform:
            #kdrew: add pseudo-count of 0.1
            i_diffrac_normalized =  np.log10(feat_df.ix[id1]['diffrac_normalized']+0.1)
        else:
            i_diffrac_normalized =  feat_df.ix[id1]['diffrac_normalized']

        zscore = (i_diffrac_normalized - mean_tmp)/std_tmp
        sliding_zscore_dict[id1] = zscore

    df = pd.DataFrame(sliding_zscore_dict.items(), columns=['ACC', 'sliding_zscore']).set_index('ACC')
    print df
    return df


def calc_fdr_correct(feat_df, unannotated_only=False):
    df = feat_df.fillna(0.0)
    if unannotated_only:
        df = df[~df.annotated]
    #df['pvalues'] = stats.norm.sf(abs(df['zscore'].values))
    df['pvalues'] = stats.norm.sf(df['zscore'].values)
    df['pvalues_fdrcor'] = fdrcorrection0(df['pvalues'])[1]

    return df[['pvalues','pvalues_fdrcor']]


def calc_sliding_fdr_correct(feat_df, unannotated_only=False):
    df = feat_df.fillna(0.0)
    if unannotated_only:
        df = df[~df.annotated]
    #df['sliding_pvalues'] = stats.norm.sf(abs(df['sliding_zscore'].values))
    df['sliding_pvalues'] = stats.norm.sf(df['sliding_zscore'].values)
    df['sliding_pvalues_fdrcor'] = fdrcorrection0(df['sliding_pvalues'])[1]

    return df[['sliding_pvalues','sliding_pvalues_fdrcor']]


if __name__ == "__main__":
    main()

