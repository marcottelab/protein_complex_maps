
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import csv
import protein_complex_maps.protein_util as pu


TINY_NUM = 0.000001
def main():

    parser = argparse.ArgumentParser(description="Plot precision recall of protein interactions for network deconvolution and correlation")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename of output plot")
    parser.add_argument("--sec_filename", action="store", dest="sec_filename", required=False, default='/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.csv',
                                            help="Filename of output plot")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
                                            help="Protein ids in which to anaylze")
    parser.add_argument("--combine", action="store_true", dest="combine", required=False,  default=False,
                                            help="Combine all profiles")
    args = parser.parse_args()

    commd_ids = ['Q86X83', 'Q9Y371', 'Q9Y6G5', 'Q9UBI1', 'Q9GZQ3', 'Q8N668',  'O60826',  'Q567U6'] 
    commd_id_map = {'Q9Y371':'sh3glb1', 'Q86X83':'commd2', 'Q9Y6G5':'commd10', 'Q9UBI1':'commd3', 'Q9GZQ3':'commd5', 'Q8N668':'commd1',  'O60826':'ccdc22',  'Q567U6':'ccdc93'} 

    dimer_ids = ['P35251', 'P12956']
    dimer_id_map = {'P12956':'P12956','P35251':'P35251'}

    #id_map = commd_id_map

    id_map = pu.get_genenames_uniprot(args.proteins)

    #kdrew: sec data (fractions 7-40) keyed on uniprot
    sec_data = dict()

    with open(args.sec_filename, 'rb') as sec_file:
        reader = csv.reader(sec_file)
        for row in reader:
            if row[0] in args.proteins:
                sec_data[ row[0] ] = map( float, row[1:34] )

    print "SEC DATA"
    print sec_data

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax3 = ax1.twinx()


    #kdrew: make initial array all equal to 1.0
    combined_data = np.repeat(1.0,len(sec_data[sec_data.keys()[0]]))
    print combined_data


    individual_color = '0.90'

    for i in args.proteins:
        #kdrew: why wouldn't it be, this is weird
        if i in args.proteins:
            if args.combine:
                ax1.plot( sec_data[i], color=individual_color, zorder=0)
                ax2.plot( sec_data[i], color=individual_color, zorder=-1)
            else:
                try:
                    ax1.plot( sec_data[i])
                    ax2.plot( sec_data[i])
                except KeyError, e:
                    print "KeyError: %s, skipping" % (e,)
                    continue

            sec_data_noised = np.array(sec_data[i]) + TINY_NUM
            combined_data = combined_data * sec_data_noised
            print "########################"
            print np.array(sec_data[i]).max() 
            print combined_data
            print combined_data.max()
            combined_data = combined_data/combined_data.max()
            print combined_data


    if args.combine:
        ax1.plot( combined_data , linewidth=4, zorder=1)
    else:
        #ax1.legend(sec_data.keys())
        ax1.legend([id_map[i] for i in args.proteins],prop={'size':7})



    xlabels = range(10,41,5)
    ax1.set_xticks(range(3,33,5))
    ax1.set_xticklabels(xlabels)
    ax1.set_xlabel('Fraction')

    toplabels = [670, 440,130,67,15]
    #ax2.set_xticks([10,13,15,19,27])
    ax2.set_xticks([9,14,16,19,26])
    ax2.set_xticklabels(toplabels, rotation='vertical')
    ax2.set_xlabel('kDa')
    ax2.set_zorder(-1)

    #ax2.set_title("Concensus SEC of Commander Complex subunits")
    #ax2.set_title("SEC of Commander Complex subunits")

    plt.savefig(args.plot_filename)



if __name__ == "__main__":
    main()




