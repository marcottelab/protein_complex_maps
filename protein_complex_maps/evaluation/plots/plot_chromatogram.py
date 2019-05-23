
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import pickle
import pandas as pd
import itertools as it

#from sklearn.metrics import precision_recall_curve
#from sklearn.metrics import average_precision_score

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'grid': False})
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

def main():

    parser = argparse.ArgumentParser(description="Plot Chromatogram")
    parser.add_argument("--chromatogram_files", action="store", dest="chromatogram_files", nargs='+', required=True, 
                                    help="Filenames of chromatograms (exported from Chromeleon)")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--labels", action="store", dest="labels", nargs='+', required=False, default=None,
                                    help="Labels for input chromatograms in order of --chromatogram_files list")
    parser.add_argument("--start_collecting", action="store", type=float, dest="start_collecting", required=False, default=None,
                                    help="Start time (min) fractions begin to be collected, (example 10.0)")
    parser.add_argument("--stop_collecting", action="store", type=float, dest="stop_collecting", required=False, default=None,
                                    help="End time (min) fractions stop being collected, (example 75.0)")
    parser.add_argument("--collection_period", action="store", type=float, dest="collection_period", required=False, default=None,
                                    help="Time period (min) fractions are collected for (example: .75 for 45 secs")
    parser.add_argument("--colors", action="store", dest="colors", nargs='+', required=False, default=["black","red","green","orange"],
                                    help="Colors to use for plotting, default=(black red green orange)")
    parser.add_argument("--ignore_lines", action="store", type=int, dest="ignore_lines", required=False, default=38,
                                    help="Number of top non-data lines to ignore in Chromeleon file. default=38")

    args = parser.parse_args()

    #kdrew: update here 04/18/19

    chromatogram_file_dict = dict()
    for i, filename in enumerate(args.chromatogram_files):

        df = pd.read_csv(filename, header=args.ignore_lines, sep='\t', thousands=',')

        if args.labels != None and len(args.labels) == len(args.chromatogram_files):
            chromatogram_file_dict[args.labels[i]] = df
        else:
            chromatogram_file_dict[filename] = df



    for i, label in enumerate(chromatogram_file_dict.keys()):
        ax = sns.lineplot(x="Time (min)", y="Value (mAU)", data=chromatogram_file_dict[label], color=it.cycle(args.colors)[i], label=label)

    #kdrew: set axis for fraction numbers
    if args.start_collecting != None or args.stop_collecting != None or args.collection_period != None:
        fraction_mins = np.arange(args.start_collecting, args.stop_collecting, args.collection_period)
        ax3 = ax.twiny()
        ax3.plot(fraction_mins, visible=False)
        ax3.set_xticks(fraction_mins[::10])
        ax3.set_xticklabels(range(len(fraction_mins))[::10])
        ax3.set_xlabel("Fractions")


    plt.legend(loc="upper right",fontsize=8)

    #plt.xlabel('Recall')
    #plt.ylabel('Precision')
    #plt.ylim([0.0, 1.05])
    #plt.xlim([0.0, 1.0])
    ##plt.title('Precision-Recall example: AUC={0:0.2f}'.format(average_precision))
    #plt.title('Precision-Recall')

    plt.savefig(args.output_file)


if __name__ == "__main__":
    main()


