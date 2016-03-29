
import argparse
import numpy as np
import operator

from scipy.stats import pearsonr

import protein_complex_maps.hhsuite.hhdomain as hhd

def main():

    parser = argparse.ArgumentParser(description="Compare complexes based on SCOP frequencies")
    parser.add_argument("--complex_scop_freq_filename", action="store", dest="scop_freq_filename", required=True, 
                                            help="Filename of scop frequency table")
    parser.add_argument("--shuffle_complex_scop_freq_filenames", nargs='+',  action="store", dest="shuffle_scop_freq_filenames", required=True, 
                                            help="Filenames of shuffled complex scop frequency tables")
    parser.add_argument("--complex_filename", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes")
    #parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
    #                                        help="Output of SCOP frequency table for each complex")
    args = parser.parse_args()

    complexes = []
    complex_file = open(args.complex_filename,"rb")

    for line in complex_file.readlines():
        complexes.append(line.split())

    complex_file.close()

    complexes_scop_freq_dict = dict()
    freq_file = open(args.scop_freq_filename,"rb")
    for line in freq_file.readlines():
        i = line.split()[0].split(':')[1]
        complexes_scop_freq_dict[i] = dict()
        for sid_freq in line.split()[1:]:
            sid = sid_freq.split(':')[0]
            freq = int(sid_freq.split(':')[1])
            complexes_scop_freq_dict[i][sid] = freq

    freq_file.close()

    shuffled_scop_freq_list = []
    for shuffle_file in args.shuffle_scop_freq_filenames:
        shuffled_complexes_scop_freq_dict = dict()
        freq_file = open(shuffle_file,"rb")
        for line in freq_file.readlines():
            i = line.split()[0].split(':')[1]
            shuffled_complexes_scop_freq_dict[i] = dict()
            for sid_freq in line.split()[1:]:
                sid = sid_freq.split(':')[0]
                freq = int(sid_freq.split(':')[1])
                shuffled_complexes_scop_freq_dict[i][sid] = freq

        freq_file.close()
        shuffled_scop_freq_list.append(shuffled_complexes_scop_freq_dict)




    complex_similarity = dict()
    complex_similarity_zscore = dict()
    #kdrew: for every complex vs every other complex, compute similarity score (dot product or correlation) of scop frequency tables
    for i in complexes_scop_freq_dict:
        complex_similarity[i] = dict()
        complex_similarity_zscore[i] = dict()
        for j in complexes_scop_freq_dict:
            similarity_metric = min_sum_score(complexes_scop_freq_dict[i], complexes_scop_freq_dict[j])
            #print "%s,%s : %s" % (i,j,similarity_metric)
            if similarity_metric > 0.0:
                shuffled_similarity_list = []
                for ii in range(len(shuffled_scop_freq_list)):
                    shuffled_similarity_list.append(min_sum_score(complexes_scop_freq_dict[i], shuffled_scop_freq_list[ii][j]))


                complex_similarity[i][j] = similarity_metric
                #kdrew: TODO: need to check for divide by zero
                complex_similarity_zscore[i][j] = (1.0*similarity_metric - np.mean(shuffled_similarity_list)) / np.std(shuffled_similarity_list)

    for i in complex_similarity_zscore:
        sorted_x = sorted(complex_similarity_zscore[i].items(), key=operator.itemgetter(1))
        #print "complex: %s" % (i,)
        for j,zscore in sorted_x:
            print "%s,%s,%s,%s,%s" % (i, j, zscore, complex_similarity[i][j], jaccard_index(complexes[int(i)],complexes[int(j)]), )

def jaccard_index(x, y):
    sx = set(x)
    sy = set(y)
    return 1.0 * len(sx.intersection(sy)) / len(sx.union(sy))


def min_sum_score(complex_freq1, complex_freq2):
    similarity_metric = 0.0
    for scop_i in complex_freq1:
        try:
            similarity_metric += min(complex_freq2[scop_i],complex_freq1[scop_i])
        except KeyError:
            continue
    for scop_j in complex_freq2:
        if scop_j not in complex_freq1.keys():
            try:
                similarity_metric += min(complex_freq2[scop_j],complex_freq1[scop_j])
            except KeyError:
                continue

    return similarity_metric
        
def dot_product(complex_freq1, complex_freq2):
    similarity_metric = 0.0
    for scop_i in complex_freq1:
        try:
            similarity_metric += complex_freq2[scop_i] * complex_freq1[scop_i]
        except KeyError:
            continue
    for scop_j in complex_freq2:
        if scop_j not in complex_freq1.keys():
            try:
                similarity_metric += complex_freq2[scop_j] * complex_freq1[scop_j]
            except KeyError:
                continue

    return similarity_metric

if __name__ == "__main__":
	main()


