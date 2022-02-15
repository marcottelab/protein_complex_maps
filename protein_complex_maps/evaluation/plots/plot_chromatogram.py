
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
#rcParams.update({'grid': False})
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
    parser.add_argument("--normalize", action="store_true", dest="normalize", required=False, default=False,
                                    help="Normalize individual chromatograms to max value = 1.0, default = False")
    parser.add_argument("--colors", action="store", dest="colors", nargs='+', required=False, default=["black","red","green","orange","pink","purple"],
                                    help="Colors to use for plotting, default=(black red green orange)")
    parser.add_argument("--standard_time", action="store", dest="standard_time", nargs='+', required=False, default=[],
                                    help="Retention time of standards run on column, example = '22.796667 31.836667 35.333333 37.683333 38.950000 39.930000 43.093333', default = ''")
    parser.add_argument("--standard_mass", action="store", dest="standard_mass", nargs='+', required=False, default=[],
                                    help="Sizes of standards run on column, example = '2000 669 443 200 150 66 29', default = ''")
    parser.add_argument("--ignore_lines", action="store", type=int, dest="ignore_lines", required=False, default=38,
                                    help="Number of top non-data lines to ignore in Chromeleon file. default=38")

    args = parser.parse_args()

    chromatogram_file_list = []
    for i, filename in enumerate(args.chromatogram_files):

        df = pd.read_csv(filename, header=args.ignore_lines, sep='\t', thousands=',')
        df['Normalized Value'] = df['Value (mAU)'] / df['Value (mAU)'].max()

        chromatogram_file_list.append(df)

    value_col_name = "Value (mAU)"
    if args.normalize:
        value_col_name = "Normalized Value"

    for i, df in enumerate(chromatogram_file_list):

        #kdrew: color is chosen based on a itertools slice of a cycle which looks a little hacky, surprised there isn't a more elegant solution
        color = [j for j in it.islice(it.cycle(args.colors),0,len(chromatogram_file_list))][i]

        label = args.chromatogram_files[i]
        if args.labels != None and len(args.labels) == len(args.chromatogram_files):
            label = args.labels[i]

        #kdrew: plot chromatogram line from dataframe
        ax = sns.lineplot(x="Time (min)", y=value_col_name, data=df, color=color, label=label)


    #kdrew: set axis for fraction numbers
    fraction_mins = None
    if args.start_collecting != None and args.stop_collecting != None and args.collection_period != None:
        fraction_mins = np.arange(args.start_collecting, args.stop_collecting, args.collection_period)
        ax3 = ax.twiny()
        #ax3.plot(fraction_mins, visible=False)
        ax3.plot(range(int(df['Time (min)'].values[0]),int(df['Time (min)'].values[-1])), visible=False)
        ax3.set_xticks(fraction_mins[::10])
        ax3.set_xticklabels(range(len(fraction_mins))[::10])
        ax3.set_xlabel("Fractions")


    std_times = [float(x) for x in args.standard_time]
    for i,x in enumerate(std_times):
        ax.axvline(x, color="black", ls='dashed', alpha=0.5, linewidth=0.5)
        ax.annotate(args.standard_mass[i], xy=(x, 1.12), xycoords='data', xytext=(-2.5, 0.0), textcoords='offset points', ha="center", rotation=90, fontsize=6,)

    print fraction_mins
    print std_times
    std_fractions = []
    for st in std_times:
        std_fractions.append(next(i for i,v in enumerate(fraction_mins) if v > st))

    print "standard_fractions"
    print ' '.join([str(x) for x in std_fractions])

    if args.normalize:
        plt.ylim(-0.1,1.18)

    plt.legend(loc="upper right",fontsize=8)

    plt.savefig(args.output_file)


if __name__ == "__main__":
    main()


