
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

    parser = argparse.ArgumentParser(description="Analyze complexes by age")
    parser.add_argument("--complex_file", action="store", dest="complex_filename", required=False, default='/home/kdrew/data/protein_complex_maps/20130828/Hs_2sp_v35.6_981cxs_complexes.txt',
                help="Filename of complex ids")
    parser.add_argument("--age_file", action="store", dest="age_filename", required=False, default='/home/kdrew/data/protein_complex_maps/20130828/tillierAnal/oma_age_2153_proteins.txt',
                help="Filename of protein age data")
    parser.add_argument("--class_file", action="store", dest="class_filename", required=False, default='/home/kdrew/data/protein_complex_maps/20130828/tillierAnal/complex_age_class.txt',
                help="Filename of protein age data")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                help="Filename of plot")
    args = parser.parse_args()


    #kdrew: read in complexes
    complexes = {}
    index = 0
    with open(args.complex_filename, 'rb') as complex_file:
        for complex_line in complex_file.readlines():
            complexes[index] = complex_line.split()
            print "index: %s proteins: %s" % (index, complexes[index])
            index+=1

    print complexes

    #kdrew: read in protein ages
    protein_age = {}
    with open(args.age_filename, 'rb') as age_file:
        for age_line in age_file.readlines():
            try:
                protein_age[age_line.split()[0]] = float(age_line.split()[1])
            except IndexError:
                print "missing value for protein age"
                continue

    print protein_age

    #kdrew: read in complex age classification (new (metazoan), mixed, old)
    complex_class = {}
    class_old = []
    class_mix = []
    class_metazoan = []
    with open(args.class_filename, 'rb') as class_file:
        for class_line in class_file.readlines():
            complex_class[class_line.split()[0]] = class_line.split()[1]
            if class_line.split()[1] == "old":
                class_old.append(int(class_line.split()[0]))
            elif class_line.split()[1] == "metazoan":
                class_metazoan.append(int(class_line.split()[0]))
            elif class_line.split()[1] == "mix":
                class_mix.append(int(class_line.split()[0]))
            else:
                print "Something wrong with class definitions"

    print complex_class
    print class_old
    print class_mix
    print class_metazoan

    complex_ages = dict()
    old_ages = []
    metazoan_ages = []
    mix_ages = []
    for cplx in complexes:
        print cplx
        complex_ages[cplx] = []
        for prot in complexes[cplx]:
            try:
                complex_ages[cplx].append(protein_age[prot])
            except KeyError:
                continue

    print complex_ages

    for cplx in class_old:
        for prot in complexes[cplx]:
            try:
                old_ages.append(protein_age[prot])
            except KeyError:
                continue

    for cplx in class_metazoan:
        for prot in complexes[cplx]:
            try:
                metazoan_ages.append(protein_age[prot])
            except KeyError:
                continue

    for cplx in class_mix:
        for prot in complexes[cplx]:
            try:
                mix_ages.append(protein_age[prot])
            except KeyError:
                continue


    print old_ages
    print mix_ages
    print metazoan_ages

    print np.mean(old_ages)
    print np.mean(mix_ages)
    print np.mean(metazoan_ages)

    #kdrew: calculate average age for each complex
    avg_age_array = []
    for cplx in complexes:
        avg_age_array.append(np.mean(complex_ages[cplx]))

    avg_age_indices = np.argsort(avg_age_array)

    #kdrew: calculate histograms for each complex
    age_hist_array = []
    for cplx in complexes:
        a = np.histogram(complex_ages[cplx], range=(0.0,0.6))[0]
        a = 1.0*a/np.sum(a)
        age_hist_array.append(a)

    age_hist_array = np.array(age_hist_array)
    #for a in age_hist_array:
    #    print a


    #kdrew: create set of ids in all complexes
    #complex_ids = set([item for sublist in complexes for item in sublist])

    

    fig = plt.figure(num=None, figsize=(8,10))
    plt.imshow(age_hist_array[avg_age_indices], aspect='auto')
    fig.savefig(args.plot_filename+"hist.pdf")

    fig = plt.figure(num=None, figsize=(8,10))
    #plt.imshow(complex_sec_average_array_woNone[mw_sort_idx_woNone], aspect='auto')
    plt.subplot(311)
    plt.hist(old_ages, range=(0.0,0.6))
    #plt.hist(old_ages, bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6])
    plt.xlim([0.0,0.6])
    plt.xlabel('Old Protein Age')
    plt.ylabel('Counts')
    plt.subplot(312)
    plt.hist(mix_ages, range=(0.0,0.6))
    #plt.hist(mix_ages, bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6]) 
    plt.xlim([0.0,0.6])
    plt.xlabel('Mixed Protein Age')
    plt.ylabel('Counts')
    plt.subplot(313)
    plt.hist(metazoan_ages, range=(0.0,0.6))
    #plt.hist(metazoan_ages, bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6])  
    plt.xlim([0.0,0.6])
    plt.xlabel('Metazoan Protein Age')
    plt.ylabel('Counts')
    fig.savefig(args.plot_filename)
    





if __name__ == "__main__":
    main()

