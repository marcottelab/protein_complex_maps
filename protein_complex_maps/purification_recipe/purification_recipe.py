
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

	args = parser.parse_args()

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

            df_dict[ msds_filename ] = df

            #for prot in args.proteins:
            #    print df[prot]



        #kdrew: score each by summing and print out top fraction in each fractionation experiment
        for exp in df_dict:
            df = df_dict[exp]

            score_mask = [ 1 if i in args.proteins else -1 for i in df.columns]

            df_scored = df * score_mask

            top_score_fraction = df_scored.sum(axis=1).idxmax()

            #print df_scored
            print "%s %s" % (top_score_fraction, df_scored.loc[top_score_fraction].sum())

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


        for x in it.product(list(df_dict[df_dict.keys()[0]].index), list(df_dict[df_dict.keys()[1]].index)):
            df1 = df_dict[df_dict.keys()[0]]
            df2 = df_dict[df_dict.keys()[1]]
            print x
            print df1.loc[x[0]]
            print df2.loc[x[1]]

            combined_df = df1.loc[x[0]].mul( df2.loc[x[1]], fill_value=0.0 )

            #kdrew: this scoring function gives a 1*abundance to proteins within in the complex and a -1*abundance to all other proteins
            #kdrew: the sum is the total score 
            score_mask_1neg1 = [ 1 if i in args.proteins else -1 for i in combined_df.index]
            df_scored = combined_df * score_mask_1neg1
            print "score %s: %s" % (x, df_scored.sum())


            #kdrew: this scoring function gives what percent the final solution will have of proteins in the complex
            score_mask_1_0 = [ 1 if i in args.proteins else 0 for i in combined_df.index]
            df_scored = combined_df * score_mask_1_0
            print "percent score %s: %s" % (x, df_scored.sum()/combined_df.sum())


if __name__ == "__main__":
	main()


