
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import multiprocessing as mp
import argparse
import cPickle
import pickle

import itertools as it

import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr


TINY_NUM = 0.00001

def main():
    parser = argparse.ArgumentParser(description="Plot abundance profile of proteins")
    parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True,
                                help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename of output plot")

    args = parser.parse_args()

    df_dict = dict()
    for msds_filename in args.msds_filenames:
        msds = pickle.load( open( msds_filename, "rb" ) )
        df_dict[msds_filename] = msds.get_data_frame()


    plot_comparison(df_dict, savefilename=args.plot_filename)

def plot_comparison(df_dict, savefilename=None, logtransform = True ):

    for (i,j) in it.combinations(df_dict.keys(),2):
        print df_dict[i]
        print df_dict[j]

        #kdrew: create union of protein ids and reindex data frames to have all protein ids and fill NaNs with 0.0
        #kdrew: effectively adds missing proteins with 0.0 values
        df_i = df_dict[i].transpose().reindex( set(df_dict[i].columns) | set(df_dict[j].columns) ).transpose().fillna(0.0)
        df_j = df_dict[j].transpose().reindex( set(df_dict[i].columns) | set(df_dict[j].columns) ).transpose().fillna(0.0)

        df_i = df_i + TINY_NUM
        df_j = df_j + TINY_NUM

        print df_i
        print df_j

        fig, ax = plt.subplots(subplot_kw=dict(axisbg='#FFFFFF'))
        ax.grid(color='#EEEEEE', linestyle='solid')
        ax.set_title("Experiment 1 vs Experiment 2", size=15)
        
        ax.set_xlabel("Experiment 1", size=15)
        ax.set_ylabel("Experiment 2", size=15)
        
        exp1_values_list = []
        exp2_values_list = []


        #kdrew: calculate cross-correlation of dataframes to find best alignment


        #kdrew: takes the min length of the array, ignores the last columns in the longer array that does not overlap
        min_len = min(len(df_j.index), len(df_i.index))
        len_diff = abs(len(df_j.index) - len(df_i.index))
        for ii in range(min_len):
            exp1_values = df_i.ix[ii] 
            exp2_values = df_j.ix[ii] 

            if logtransform:
                exp1_values = np.log( exp1_values )
                exp2_values = np.log( exp2_values ) 

            scatter = ax.scatter( exp1_values, exp2_values ) 
            exp1_values_list.append(exp1_values)
            exp2_values_list.append(exp2_values)


        fig.savefig(savefilename)

        print pd.concat(exp1_values_list)
        print pd.concat(exp2_values_list)

        print "pearsonr: %s: %s" % (pearsonr(pd.concat(exp1_values_list), pd.concat(exp2_values_list)))
        print "spearmanr: %s: %s" % (spearmanr(pd.concat(exp1_values_list), pd.concat(exp2_values_list)))
        

if __name__ == "__main__":
	main()



