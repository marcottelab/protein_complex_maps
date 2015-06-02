
import logging
import numpy as np
import os.path
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy
import pylab
import itertools as it

import numpy.linalg as linalg

from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

import sklearn.metrics as skm

import protein_complex_maps.normalization_util as nu
import protein_complex_maps.read_data as rd
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.plots.plot_profile as pp
import protein_complex_maps.bicluster_generator as bg
import protein_complex_maps.annealer as anl
import protein_complex_maps.physical_interaction_clustering as pic
import protein_complex_maps.external.npeet.entropy_estimators as ee


import argparse
import pickle
import scipy.cluster.hierarchy as sch

import protein_complex_maps.protein_util as pu

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

	parser = argparse.ArgumentParser(description="Create purification recipe")
	parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True, 
                                    help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
                                    help="Protein ids in which to anaylze")
	parser.add_argument("--protein_percent_threshold", action="store", dest="protein_percent_threshold", type=float, required=False, default=0.9,
                                    help="Percentage of proteins that need to be present in input fractions (expressed as decimal), default 0.9 (ie 90%)")

	args = parser.parse_args()

        results_list = []

        df_dict = dict()
        #kdrew: store data in data frames
        for msds_filename in args.msds_filenames:
            msds = pickle.load( open( msds_filename, "rb" ) )
            df = msds.get_data_frame()

            print "filename: %s" % msds_filename

            missing_proteins = [ i for i in args.proteins if i not in df.columns ]
            print "missing: %s" % missing_proteins
            #kdrew: add in 0.0 rows for missing entries
            for missing_prot in missing_proteins:
                df[missing_prot] = np.repeat( 0.0, len(df.index) )


            ##kdrew: normalize by fractions (max of all fraction vectors = 1.0)
            #df = df.div( df.max(axis=1), axis=0 )

            #kdrew: normalize by fractions (sum of all fraction vectors = 1.0, i.e. probability)
            df = df.div( df.sum(axis=1), axis=0 )

            #kdrew: count how many input proteins are present (> 0.0)
            #protein_count = (combined_df > 0.0).sum()
            #print "protein_count %s / %s = %s" % (protein_count, len(args.proteins), (1.0*protein_count)/len(args.proteins))

            #kdrew: threshold out fraction vectors that do not have enough of the input proteins
            percent_vector = 1.0*(df[args.proteins] > 0.0).sum(axis=1)/len(args.proteins)
            df = df[percent_vector > args.protein_percent_threshold]
            print df

            df_dict[ msds_filename ] = df

            #for prot in args.proteins:
            #    print df[prot]



        #kdrew: testing
        #kdrew: score each by summing and print out top fraction in each fractionation experiment
        for exp in df_dict:
            df = df_dict[exp]

            score_mask = [ 1 if i in args.proteins else -1 for i in df.columns]

            df_scored = df * score_mask

            try:
                top_score_fraction = df_scored.sum(axis=1).idxmax()

                #print df_scored
                print "%s %s" % (top_score_fraction, df_scored.loc[top_score_fraction].sum())
            except ValueError:
                print "no fraction passed threshold"

        #kdrew: testing single experiment
        for x in list(df_dict[df_dict.keys()[0]].index):
            df1 = df_dict[df_dict.keys()[0]]
            print "df1: %s" % df1
            print "x: %s" % x
            print "df1.loc[x] %s" % df1.loc[x]

            for prot in args.proteins:
                print "%s : %s" % (prot, df1.loc[x][prot])

            #kdrew: this scoring function gives a 1*abundance to proteins within in the complex and a -1*abundance to all other proteins
            #kdrew: the sum is the total score 
            score_mask_1neg1 = [ 1 if i in args.proteins else -1 for i in df1.loc[x].index]
            df_scored = df1.loc[x] * score_mask_1neg1
            print "score %s: %s" % (x, df_scored.sum())


            #kdrew: this scoring function gives what percent the final solution will have of proteins in the complex
            score_mask_1_0 = [ 1 if i in args.proteins else 0 for i in df1.loc[x].index]
            print score_mask_1_0
            df_scored = df1.loc[x] * score_mask_1_0
            print "sum %s: %s %s" % (x, df_scored.sum(), df1.loc[x].sum())
            print "percent score %s: %s" % (x, (df_scored.sum()/df1.loc[x].sum()))
            #print "df_scored: %s" % df_scored
            #print "percent score %s: %s" % (x, (df_scored.sum()))

        #kdrew: testing another single experiment
        for x in it.product(list(df_dict[df_dict.keys()[1]].index)):
            df1 = df_dict[df_dict.keys()[1]]
            print "df1: %s" % df1
            print "x: %s" % x
            print "df1.loc[x[0]] %s" % df1.loc[x[0]]

            for prot in args.proteins:
                print "%s : %s" % (prot, df1.loc[x[0]][prot])

            #kdrew: this scoring function gives a 1*abundance to proteins within in the complex and a -1*abundance to all other proteins
            #kdrew: the sum is the total score 
            score_mask_1neg1 = [ 1 if i in args.proteins else -1 for i in df1.loc[x[0]].index]
            df_scored = df1.loc[x[0]] * score_mask_1neg1
            print "score %s: %s" % (x[0], df_scored.sum())


            #kdrew: this scoring function gives what percent the final solution will have of proteins in the complex
            score_mask_1_0 = [ 1 if i in args.proteins else 0 for i in df1.loc[x].index]
            print score_mask_1_0
            df_scored = df1.loc[x] * score_mask_1_0
            print "sum %s: %s %s" % (x, df_scored.sum(), df1.loc[x].sum())
            print "percent score %s: %s" % (x, (df_scored.sum()/df1.loc[x].sum()))
            #print "df_scored: %s" % df_scored
            #print "percent score %s: %s" % (x, (df_scored.sum()))





        for r in range(1, len(df_dict)+1):
            for exp_list in it.combinations( df_dict.keys(), r ):

                df_list = [ df_dict[exp] for exp in exp_list ] 
                df_index_list = [ list(df.index) for df in df_list ]
                print "combinations: %s" % (','.join(exp_list))
                print "possible fractions: %s" % (df_list,)

                for fractions in it.product(*df_index_list):

                    print fractions
                    #print df_list[0].loc[fractions[0]]
                    #print df_list[1].loc[fractions[1]]

                    #kdrew: initialize vector to be the first fraction
                    combined_vector = df_list[0].loc[fractions[0]]
                    #kdrew: multiply additional fraction vectors
                    for i in range(1,len(fractions)):
                        combined_vector = combined_vector .mul( df_list[i].loc[fractions[i]], fill_value=0.0 )

                    pr = score_fraction(combined_vector, args.proteins)

                    pr.experiments=exp_list
                    pr.fractions=fractions

                    results_list.append(pr)

        for result in results_list:
            print result


def score_fraction(vector, proteins):

    #kdrew: this scoring function gives a 1*abundance to proteins within in the complex and a -1*abundance to all other proteins
    #kdrew: the sum is the total score 
    score_mask_1neg1 = [ 1 if i in proteins else -1 for i in vector.index]
    vector_scored = vector * score_mask_1neg1
    print "score %s" % (vector_scored.sum())


    #kdrew: count how many original proteins are still present (> 0.0)
    protein_count = (vector[proteins] > 0.0).sum()
    proteins_present = vector.index[(vector[proteins] > 0.0)].tolist()
    protein_percent = (1.0*protein_count)/len(proteins)
    print "protein_count %s / %s = %s" % (protein_count, len(proteins), protein_percent)

    #kdrew: this scoring function gives what percent the final solution will have of proteins in the complex
    score_mask_1_0 = [ 1 if i in proteins else 0 for i in vector.index]
    vector_scored = vector * score_mask_1_0
    purity_percent = vector_scored.sum()/vector.sum()
    print "percent score %s" % (purity_percent)

    pr = PurificationRecipe(proteins=proteins_present, protein_percent=protein_percent, purity_percent=purity_percent)

    return pr


class PurificationRecipe(object):

    def __init__(self, experiments=[], fractions=[], proteins=[], protein_percent=None, purity_percent=None):
        self.experiments = experiments
        self.fractions = fractions
        self.proteins = proteins
        self.protein_percent = protein_percent
        self.purity_percent = purity_percent

    def __str__(self,):
        return "fractions: %s, protein_percent: %s, purity_percent: %s" % (self.fractions, self.protein_percent, self.purity_percent)

if __name__ == "__main__":
	main()


