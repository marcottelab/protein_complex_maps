
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats

from sklearn.ensemble import RandomForestClassifier

from pyemd import emd

import protein_complex_maps.features.ExtractFeatures.Features as eff

def main():

    parser = argparse.ArgumentParser(description="Calculate difference features between two fractionation experiments")
    parser.add_argument("--elution_files", action="store", nargs='+', dest="elution_files", required=True, 
                                    help="Elution files (.elut)")
    parser.add_argument("--features", action="store", nargs='+', dest="features", required=False, default=['diffrac'],
                                    help="Features to calculate: diffrac (L1-norm of difference), diffrac_percent, diffrac_normalized, pearsonr, poisson, mean_abundance, emd")
    parser.add_argument("--classifiers", action="store", nargs='+', dest="classifiers", required=False, default=[],
                                    help="Classifiers to train: random_forest")
    parser.add_argument("--training_labels", action="store", dest="training_labels", required=False, default=None, 
                                    help="Filename of training labels, default=None")
    parser.add_argument("--training_id_column", action="store", dest="training_id_column", required=False, default='GENE PRODUCT ID', 
                                    help="Column name of id in training label file, default='GENE PRODUCT ID'")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")

    args = parser.parse_args()
    
    elutions = []
    for efile in args.elution_files:
        elut = eff.Elut()
        elut.load(efile,format='tsv')
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

        #kdrew: read in training labels for machine learning scores
        if args.training_labels != None:
            #kdrew: add in training labels
            training_labels_df = pd.read_table(args.training_labels)
            labels = [i in training_labels_df[args.training_id_column].values for i in feature_df.index]
            feature_df['label'] = labels

            if 'random_forest' in args.classifiers:
                return_df = random_forest(feature_df, args.features)
                feature_df = join_feature(feature_df,return_df)

        if args.out_filename != None:
            feature_df.sort_values(args.features[0]).to_csv(args.out_filename)
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


def random_forest(feature_df, features, random_state=0, n_estimators=1500):

    feature_nonan_df = feature_df.fillna(0.0)
    #kdrew: create training set by randomly sampling half of the proteins
    feature_training_df = feature_nonan_df.sample(len(feature_nonan_df)/2, random_state=random_state)

    #kdrew: get all of the entries that are not in the training set
    feature_test_df = feature_nonan_df.loc[[x not in feature_training_df.index for x in feature_nonan_df.index]][:]

    feature_training_data = feature_training_df[features].values
    feature_training_labels = feature_training_df['label'].values
    feature_test_data = feature_test_df[features].values
    feature_test_labels = feature_test_df['label'].values

    RF_model =  RandomForestClassifier(random_state=random_state, n_estimators=n_estimators, oob_score=True) #, class_weight="balanced")
    RF_model.fit(feature_training_data, feature_training_labels)
    test_rf_preds = RF_model.predict_proba(feature_test_data)
    training_rf_preds = RF_model.predict_proba(feature_training_data)

    print list(RF_model.classes_)
    true_index = list(RF_model.classes_).index(True)
    feature_test_df['rf_true_prob'] = [x[true_index] for x in test_rf_preds]
    feature_test_df['rf_test_training'] = 'test'
    feature_training_df['rf_true_prob'] = [x[true_index] for x in training_rf_preds]
    feature_training_df['rf_test_training'] = 'train'

    feature_test_training_df = pd.concat([feature_test_df,feature_training_df])
    return feature_test_training_df[['rf_test_training','rf_true_prob']]



if __name__ == "__main__":
    main()

