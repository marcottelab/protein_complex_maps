
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

def main():

    parser = argparse.ArgumentParser(description="Analyze complexes with SEC data")
    parser.add_argument("--complex_file", action="store", dest="complex_filename", required=False, default='/home/kdrew/data/corum/allComplexesCore_acc.txt',
                help="Filename of complex ids")
    parser.add_argument("--sec_file", action="store", dest="sec_filename", required=False, default='/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.csv',
                help="Filename of SEC data")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                help="Filename of plot")
    parser.add_argument("--uniprot_seq_pickle", action="store", dest="seq_pickle", required=False, default=None,
                help="Filename of uniprot seq pickle")
    args = parser.parse_args()


    SEC_START_COL = 1
    SEC_END_COL = 35

    #kdrew: read in complexes
    complexes = []
    with open(args.complex_filename, 'rb') as complex_file:
        for complex_line in complex_file.readlines():
            complexes.append(complex_line.split())

    print complexes

    #kdrew: create set of ids in all complexes
    complex_ids = set([item for sublist in complexes for item in sublist])

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

    if args.seq_pickle == None:
        complex_seq_map = pu.get_sequences_uniprot(list(complex_ids))
        pickle.dump(complex_seq_map, open('./corum_uniprot_seq.p','wb'))

    else:
        complex_seq_map = pickle.load(open(args.seq_pickle, 'rb'))
    
    print complex_seq_map

    
    #kdrew: calculate molecular weight for total complex, index is same for complexes and mw list
    complex_mw = []
    for cmplx in complexes:
        weight = 0.0
        for c_id in cmplx:
            try:
                weight += su.molecular_weight(complex_seq_map[c_id], seq_type="protein")
            except (KeyError, TypeError) :
                print "Missing entry, ignore complex"
                weight = None
                break

        print "%s : %s" % (cmplx, weight)
        complex_mw.append(weight)


    #kdrew: calculate average SEC profile for complex, index same as complexes list
    complex_sec_average = []
    complex_sec_weight_fraction = []
    complex_sec_max_fraction = []
    for cmplx in complexes:
        sec_cnt = 0.0
        complex_sec_array = np.zeros(SEC_END_COL - SEC_START_COL)
        for c_id in cmplx:
            try:
                complex_sec_array = complex_sec_array + np.array(sec_data[c_id])
                sec_cnt+=1
            except KeyError:
                continue

        complex_sec_array = 1.0*complex_sec_array/sec_cnt
        complex_sec_average.append(complex_sec_array)

        #kdrew: average fraction index weighted by value in sec 
        complex_sec_weight_fraction.append(np.average(range(SEC_START_COL,SEC_END_COL), weights=complex_sec_array))   
        complex_sec_max_fraction.append(np.argmax(complex_sec_array))   

    print complex_sec_average
    complex_sec_average_array = np.array(complex_sec_average)

    complex_mw_array = np.array(complex_mw)
    complex_sec_weight_fraction_array = np.array(complex_sec_weight_fraction)
    complex_sec_max_fraction_array = np.array(complex_sec_max_fraction)
    mw_sort_idx = np.argsort(complex_mw_array)

    print complex_sec_average_array[mw_sort_idx]
    print complex_mw_array[mw_sort_idx]
    print complex_sec_weight_fraction_array[mw_sort_idx]

    nonNone_index = [i for i in range(len(complex_mw_array)) if complex_mw_array[i] != None and ~np.isnan(complex_sec_weight_fraction_array[i]) ]
    print nonNone_index

    print s.spearmanr(complex_mw_array, complex_sec_weight_fraction_array)
    print s.spearmanr(complex_mw_array[nonNone_index], complex_sec_weight_fraction_array[nonNone_index])

    print s.spearmanr(complex_mw_array, complex_sec_max_fraction_array)
    print s.spearmanr(complex_mw_array[nonNone_index], complex_sec_max_fraction_array[nonNone_index])

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
    






if __name__ == "__main__":
    main()

