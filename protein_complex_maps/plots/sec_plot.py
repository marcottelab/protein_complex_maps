
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import sklearn.mixture
from scipy.optimize import curve_fit

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
    parser.add_argument("--average", action="store_true", dest="average", required=False,  default=False,
                                            help="Average all profiles")
    args = parser.parse_args()

    #commd_ids = ['Q86X83', 'Q9Y371', 'Q9Y6G5', 'Q9UBI1', 'Q9GZQ3', 'Q8N668',  'O60826',  'Q567U6'] 
    #commd_id_map = {'Q9Y371':'sh3glb1', 'Q86X83':'commd2', 'Q9Y6G5':'commd10', 'Q9UBI1':'commd3', 'Q9GZQ3':'commd5', 'Q8N668':'commd1',  'O60826':'ccdc22',  'Q567U6':'ccdc93'} 

    #dimer_ids = ['P35251', 'P12956']
    #dimer_id_map = {'P12956':'P12956','P35251':'P35251'}

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
    #ax2 = ax1.twiny()
    #ax3 = ax1.twinx()


    #kdrew: make initial array all equal to 1.0
    combined_data = np.repeat(1.0,len(sec_data[sec_data.keys()[0]]))
    print combined_data
    total_data = np.repeat(0.0,len(sec_data[sec_data.keys()[0]]))


    individual_color = '0.90'

    proteins_present = []
    for i in args.proteins:
        #kdrew: why wouldn't it be, this is weird
        if i in args.proteins:
            if args.combine:
                ax1.plot( sec_data[i], color=individual_color, zorder=0)
                #ax2.plot( sec_data[i], color=individual_color, zorder=-1)
            elif args.average:
                ax1.plot( sec_data[i], color=individual_color, zorder=0)
            else:
                try:
                    ax1.plot( sec_data[i])
                    #ax2.plot( sec_data[i])
                    proteins_present.append(i)
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
            total_data = total_data + sec_data[i]
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print total_data
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"



    if args.combine:
        ax1.plot( combined_data , linewidth=4, zorder=1)
    if args.average:
        ax1.plot( total_data/len(args.proteins), linewidth=4, zorder=1)
    else:
        #ax1.legend(sec_data.keys())
        ax1.legend([id_map[i] for i in proteins_present],prop={'size':7})
        print [id_map[i] for i in proteins_present]
        print id_map
        #print "not setting legend"



    #xlabels = range(10,41,5)
    #ax1.set_xticks(range(3,33,5))
    #ax1.set_xticks([3.0,8.0,13.0,18.0,23.0,28.0])
    #ax1.set_xticklabels(xlabels)
    xlabels = ['','','','10','','','','','15','','','','','20','','','','','25','','','','','30','','','','','35','','','','']
    ax1.set_xticks(range(0,33))
    ax1.set_xticklabels(xlabels)
    ax1.set_xlabel('Fraction')

    toplabels = ['','','','','','','','','','670','','','','','440','','130','','','67','','','','','','','15','','','','','','']
    #toplabels = [670,440,130,67,15]
    #ax2.set_xticks([10,13,15,19,27])
    #ax2.set_xticks([9,14,16,19,26])
    #ax2.set_xticks(range(0,33))
    #ax2.set_xticklabels(toplabels, rotation='vertical')
    #ax2.set_xlabel('kDa')
    #ax2.set_zorder(-1)

    #ax2.set_title("Concensus SEC of Commander Complex subunits")
    #ax2.set_title("SEC of Commander Complex subunits")

    plt.savefig(args.plot_filename)


    print total_data
    p0 = [1., 12.0, 1.]
    bins = np.arange(0,len(sec_data[sec_data.keys()[0]]))
    print bins
    #kdrew: mask is for commander complex
    mask_total_data =  total_data * [0.0,0.0,0.0,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.0]
    coeff, var_matrix = curve_fit(gauss, bins, mask_total_data, p0=p0)

    # Get the fitted curve
    #hist_fit = gauss(bin_centres, *coeff)

    #plt.plot(bin_centres, hist, label='Test data')
    #plt.plot(bin_centres, hist_fit, label='Fitted data')

    # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
    print 'masked Raw Fitted mean = ', coeff[1]
    print 'masked Adjusted Fitted mean = ', coeff[1] + 7
    print 'masked Fitted standard deviation = ', coeff[2]

#http://stackoverflow.com/questions/11507028/fit-a-gaussian-function
# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

if __name__ == "__main__":
    main()




