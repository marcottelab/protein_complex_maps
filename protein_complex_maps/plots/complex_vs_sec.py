
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import pickle
import csv


import scipy.stats as s
import Bio.SeqUtils as su

import protein_complex_maps.protein_util as pu

TINY_NUM = 0.000001

def main():

    parser = argparse.ArgumentParser(description="Analyze complexes with SEC data")
    parser.add_argument("--complex_file", action="store", dest="complex_filename", required=False, default='/home/kdrew/data/protein_complex_maps/20130828/Hs_2sp_v35.6_981cxs_complexes.txt',
                help="Filename of complex ids")
    parser.add_argument("--sec_file", action="store", dest="sec_filename", required=False, default='/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.csv',
                help="Filename of SEC data")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                help="Filename of plot")
    parser.add_argument("--uniprot_map_pickle", action="store", dest="map_pickle", required=False, default=None,
                help="Filename of uniprot map pickle")
    parser.add_argument("--uniprot_seq_pickle", action="store", dest="seq_pickle", required=False, default=None,
                help="Filename of uniprot seq pickle")
    parser.add_argument("--nomap", action="store_true", dest="nomap", required=False, default=False,
                help="Do not map input ids from ENSG to ACC, corum file already has ACC")
    args = parser.parse_args()


    SEC_START_COL = 1
    SEC_END_COL = 35

    #sec_filename = '/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.csv'

    #subcomplex_ids = dimer_ids
    #id_map = dimer_id_map
    #plot_filename = "complex_vs_sec.pdf"
    #combine=False

    
    #kdrew: read in complexes
    complexes = []
    with open(args.complex_filename, 'rb') as complex_file:
        for complex_line in complex_file.readlines():
            complexes.append(complex_line.split())

    print complexes

    #kdrew: create set of ids in all complexes
    complex_ids = set([item for sublist in complexes for item in sublist])

    #kdrew: map ENSEMBL ids from complex file to uniprot ids (SEC data has uniprot ids)
    #kdrew: pickle results so not to hammer uniprot's webservice
    if args.map_pickle == None:
        #kdrew: map ids from ENSEMBL_ID to ACC
        complex_id_uniprot_map = pu.map_protein_ids(complex_ids, 'ENSEMBL_ID', 'ACC', reviewed=True)

        pickle.dump(complex_id_uniprot_map, open('./complex_uniprot_map.p','wb'))

    else:
        complex_id_uniprot_map = pickle.load(open(args.map_pickle, 'rb'))
    

    #kdrew: flatten dictionary by taking for first entry in list to be value
    complex_id_uniprot_single_map = {key : complex_id_uniprot_map[key][0] for key in complex_id_uniprot_map if len(complex_id_uniprot_map[key]) > 0}
    print complex_id_uniprot_map

    #kdrew: sec data (fractions 7-40) keyed on uniprot
    sec_data = dict()

    #kdrew: readin sec data
    header_count = 0
    with open(args.sec_filename, 'rb') as sec_file:
        reader = csv.reader(sec_file)
        for row in reader:
            #kdrew: file has two header lines
            if header_count < 2:
                header_count+=1
            else:
                #kdrew: uniprot id is in entry 0, data is in entries 1-35, remaining data is stdev + misc
                sec_data[ row[0] ] = map( float, row[SEC_START_COL:SEC_END_COL] )

    #print sec_data

    #kdrew: map sec data to complex ids
    complex_sec_data = dict()
    complex_sec_ids = dict()
    for c_id in complex_ids:
        print c_id 
        for uni_id in complex_id_uniprot_map[c_id]:
            try:
                complex_sec_data[c_id] = sec_data[uni_id]
                complex_sec_ids[c_id] = uni_id
                break
            except KeyError:
                continue

    if args.seq_pickle == None:
        complex_seq_map = pu.get_sequences_uniprot(complex_id_uniprot_single_map.values())
        pickle.dump(complex_seq_map, open('./complex_uniprot_seq.p','wb'))

    else:
        complex_seq_map = pickle.load(open(args.seq_pickle, 'rb'))
    
    print complex_seq_map

    
    #for seq_id in complex_seq_map:
    #    print complex_seq_map[seq_id]
    #    print "%s : %s" % (seq_id, su.molecular_weight(complex_seq_map[seq_id], seq_type="protein"))

    print len(complex_sec_data)
    print len(complex_ids)
    print len(complex_seq_map)

    #kdrew: calculate molecular weight for total complex, index is same for complexes and mw list
    complex_mw = []
    for cmplx in complexes:
        weight = 0.0
        for c_id in cmplx:
            try:
                weight += su.molecular_weight(complex_seq_map[complex_id_uniprot_single_map[c_id]], seq_type="protein")
            except KeyError:
                print "Missing entry, ignore complex"
                weight = None
                break

        print "%s : %s" % (cmplx, weight)
        complex_mw.append(weight)


    #kdrew: calculate average SEC profile for complex, index same as complexes list
    complex_sec_average = []
    complex_sec_weight_fraction = []
    for cmplx in complexes:
        sec_cnt = 0.0
        complex_sec_array = np.zeros(SEC_END_COL - SEC_START_COL)
        for c_id in cmplx:
            try:
                complex_sec_array = complex_sec_array + np.array(complex_sec_data[c_id])
                sec_cnt+=1
            except KeyError:
                continue

        complex_sec_array = 1.0*complex_sec_array/sec_cnt
        complex_sec_average.append(complex_sec_array)

        #kdrew: average fraction index weighted by value in sec 
        complex_sec_weight_fraction.append(np.average(range(SEC_START_COL,SEC_END_COL), weights=complex_sec_array))   

    print complex_sec_average
    complex_sec_average_array = np.array(complex_sec_average)

    complex_mw_array = np.array(complex_mw)
    complex_sec_weight_fraction_array = np.array(complex_sec_weight_fraction)
    mw_sort_idx = np.argsort(complex_mw_array)

    print complex_sec_average_array[mw_sort_idx]
    print complex_mw_array[mw_sort_idx]
    print complex_sec_weight_fraction_array[mw_sort_idx]

    nonNone_index = [i for i in range(len(complex_mw_array)) if complex_mw_array[i] != None and ~np.isnan(complex_sec_weight_fraction_array[i]) ]
    print nonNone_index

    print s.spearmanr(complex_mw_array, complex_sec_weight_fraction_array)
    print s.spearmanr(complex_mw_array[nonNone_index], complex_sec_weight_fraction_array[nonNone_index])

    complex_mw_array_woNone = complex_mw_array[nonNone_index]
    complex_sec_average_array_woNone = complex_sec_average_array[nonNone_index]

    mw_sort_idx_woNone = np.argsort(complex_mw_array_woNone)

    fig = plt.figure(num=None, figsize=(8,10))
    #plt.matshow(complex_sec_average_array[mw_sort_idx])
    plt.imshow(complex_sec_average_array_woNone[mw_sort_idx_woNone], aspect='auto')
    plt.xlabel('SEC Fraction')
    plt.ylabel('Complex rank (ordered by increasing molecular weight)')
    plt.colorbar()
    fig.savefig(args.plot_filename)
    




def plot_profile():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax3 = ax1.twinx()


    #kdrew: make initial array all equal to 1.0
    combined_data = np.repeat(1.0,len(sec_data[sec_data.keys()[0]]))
    print combined_data


    individual_color = '0.90'

    #plt.plot( sec_data['Q16401'] )
    #for i in sec_data.keys():
    for i in subcomplex_ids:
            if i in subcomplex_ids:
                    if combine:
                            ax1.plot( sec_data[i], color=individual_color, zorder=0)
                            ax2.plot( sec_data[i], color=individual_color, zorder=-1)
                    else:
                            ax1.plot( sec_data[i])
                            ax2.plot( sec_data[i])

                    sec_data_noised = np.array(sec_data[i]) + TINY_NUM
                    combined_data = combined_data * sec_data_noised
                    print "########################"
                    print np.array(sec_data[i]).max() 
                    print combined_data
                    print combined_data.max()
                    combined_data = combined_data/combined_data.max()
                    print combined_data


    if combine:
            ax1.plot( combined_data , linewidth=4, zorder=1)
    else:
            #ax1.legend(sec_data.keys())
            ax1.legend([id_map[i] for i in subcomplex_ids])

    xlabels = range(10,41,5)

    toplabels = [670, 440,130,67,15]

    ax1.set_xticks(range(4,33,5))
    ax1.set_xticklabels(xlabels)
    ax1.set_xlabel('Fraction')

    ax2.set_xticks([10,13,15,19,27])
    ax2.set_xticklabels(toplabels)
    ax2.set_xlabel('kDa')
    ax2.set_zorder(-1)

    #ax2.set_title("Concensus SEC of Commander Complex subunits")
    #ax2.set_title("SEC of Commander Complex subunits")

    plt.savefig(plot_filename)



if __name__ == "__main__":
    main()

