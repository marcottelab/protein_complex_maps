
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

import protein_complex_maps.preprocessing_util.complexes.complex_merge as cm

def main():

    parser = argparse.ArgumentParser(description="Split gold standard complexes into test and training")
    parser.add_argument("--input_complexes", action="store", dest="complex_filename", required=True, 
                                    help="Filename of gold standard complexes")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", required=False, default=None,
                                    help="Random seed for shuffling, default=None")
    parser.add_argument("--size_threshold", action="store", type=int, dest="size_threshold", required=False, default=None,
                                    help="Threshold for complex size (suggested 30), default=None")
    parser.add_argument("--remove_large_complexes", action="store_true", dest="remove_large_complexes", required=False, default=False,
                                    help="Remove complexes above threshold, default=False")
    parser.add_argument("--threshold_fraction", action="store", type=float, dest="threshold_fraction", required=False, default=None,
                                    help="Fraction of pairs to include from size thresholded complexes (suggested 0.1), default=None")
    parser.add_argument("--merge_threshold", action="store", type=float, dest="merge_threshold", required=False, default=1.0,
                                            help="Jiccard similarity threshold on which to merge, default=1.0 (remove exact copies)")
    parser.add_argument("--remove_largest", action="store_true", dest="remove_largest", required=False, default=False,
                                            help="Instead of merging similar clusters, remove largest")

    args = parser.parse_args()

    if args.random_seed != None:
        print "setting random seed"
        np.random.seed(args.random_seed)

    #kdrew: origin of this code is from ipython notebook: corum_test_train_revisit, heavily refactored and modified

    complexes = set()
    f = open(args.complex_filename,"rb")
    for line in f.readlines():
        #kdrew: ignore singletons
        if len(line.split()) > 1:
            complexes.add(frozenset(line.split()))

    merged_complexes = cm.merge_complexes(complexes, args.merge_threshold, remove_largest=args.remove_largest)

    test_list_clean, train_list_clean, test_ppis_clean, train_ppis_clean, neg_test_ppis_clean, neg_train_ppis_clean = split_complexes(merged_complexes, size_threshold=args.size_threshold, threshold_fraction=args.threshold_fraction, remove_large_complexes=args.remove_large_complexes, input_complexes=complexes)

    test_list_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.test.txt'
    train_list_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.train.txt' 
    test_ppis_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.test_ppis.txt'
    train_ppis_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.train_ppis.txt' 
    neg_test_ppis_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.neg_test_ppis.txt'  
    neg_train_ppis_clean_filename = '.'.join(args.complex_filename.split('.')[:-1]) + '.neg_train_ppis.txt'  

    #kdrew: write sets to outfiles
    test_list_out = open(test_list_clean_filename,"wb")
    for c in test_list_clean:
        test_list_out.write(" ".join(c)+"\n")
    test_list_out.close()

    train_list_out = open(train_list_clean_filename,"wb")
    for c in train_list_clean:
        train_list_out.write(" ".join(c)+"\n")
    train_list_out.close()

    test_out = open(test_ppis_clean_filename,"wb")
    for ppi in test_ppis_clean:
        #test_out.write("%s\t%s\n" % (list(ppi)[0], list(ppi)[1]))
        test_out.write("%s\t%s\n" % (tuple(sorted(ppi))))
    test_out.close()

    train_out = open(train_ppis_clean_filename,"wb")
    for ppi in train_ppis_clean:
        #train_out.write("%s\t%s\n" % (list(ppi)[0], list(ppi)[1]))
        train_out.write("%s\t%s\n" % (tuple(sorted(ppi))))
    train_out.close()

    neg_test_out = open(neg_test_ppis_clean_filename,"wb")
    for ppi in neg_test_ppis_clean:
        #neg_test_out.write("%s\t%s\n" % (list(ppi)[0], list(ppi)[1]))
        neg_test_out.write("%s\t%s\n" % (tuple(sorted(ppi))))
    neg_test_out.close()

    neg_train_out = open(neg_train_ppis_clean_filename,"wb")
    for ppi in neg_train_ppis_clean:
        #neg_train_out.write("%s\t%s\n" % (list(ppi)[0], list(ppi)[1]))
        neg_train_out.write("%s\t%s\n" % (tuple(sorted(ppi))))
    neg_train_out.close()


#kdrew: split complexes into test and training, return orthogonal ppi sets
def split_complexes(complexes, size_threshold=None, threshold_fraction=None, remove_large_complexes=False, input_complexes=None): 

    #kdrew: randomly shuffle complexes 
    np.random.shuffle(complexes)
    split_id = len(complexes)/2
    test_list = complexes[:split_id]
    train_list = complexes[split_id:]

    test_list_clean, train_list_clean = remove_overlapping_complexes(test_list, train_list)

    #kdrew: for complexes in test_list, generate all pairs
    test_ppis = set()
    #kdrew: holds all pairs regardless of thresholding
    full_test_ppis = set()
    for clust in test_list:
        #print "test clust"
        #print clust
        #kdrew: for large complexes sample only a fraction of pairs
        if size_threshold != None and threshold_fraction != None and len(clust) > size_threshold:
            large_pairs = [pair for pair in it.combinations( set( clust ), 2)]
            np.random.shuffle(large_pairs)
            #kdrew: calculate the index of a fraction of the list
            fraction_index = int(len(large_pairs)*threshold_fraction)
            for pair in large_pairs[:fraction_index]:
                test_ppis.add(frozenset(pair))
        else:
            for pair in it.combinations( set( clust ), 2):
                test_ppis.add(frozenset(pair))

        #kdrew: holds all pairs regardless of thresholding
        for pair in it.combinations( set( clust ), 2):
            full_test_ppis.add(frozenset(pair))
    
    #kdrew: for complexes in train_list, generate all pairs
    train_ppis = set()
    #kdrew: holds all pairs regardless of thresholding, but only of input complexes not what was merged or thrownout previously (need to fix)
    full_train_ppis = set()
    for clust in train_list:
        #print "train clust"
        #print clust
        #kdrew: for large complexes sample only a fraction of pairs
        if size_threshold != None and threshold_fraction != None and len(clust) > size_threshold:
            large_pairs = [pair for pair in it.combinations( set( clust ), 2)]
            np.random.shuffle(large_pairs)
            #kdrew: calculate the index of a fraction of the list
            fraction_index = int(len(large_pairs)*threshold_fraction)
            for pair in large_pairs[:fraction_index]:
                train_ppis.add(frozenset(pair))
        else:
            for pair in it.combinations( set( clust ), 2):
                train_ppis.add(frozenset(pair))

        #kdrew: holds all pairs regardless of thresholding
        for pair in it.combinations( set( clust ), 2):
            full_train_ppis.add(frozenset(pair))
    

    #kdrew: remove intersection from test and train ppi sets
    test_ppis_clean, train_ppis_clean = remove_intersection_random(test_ppis, train_ppis)

    #kdrew: remove test_ppis that show up in training list
    test_ppis_overlap = set()
    for clust in train_list_clean:
        #print clust
        for test_ppi in test_ppis_clean:
            #print test_ppi
            if test_ppi in [set(x) for x in it.combinations(clust,2)]:
                test_ppis_overlap.add(test_ppi)
    print "removing test_ppis overlapped with training clusters: %s" % len(test_ppis_overlap)
    test_ppis_clean = test_ppis_clean - test_ppis_overlap

    #kdrew: remove train_ppis that show up in test list
    train_ppis_overlap = set()
    for clust in test_list_clean:
        #print clust
        for train_ppi in train_ppis_clean:
            #print train_ppi
            if train_ppi in [set(x) for x in it.combinations(clust,2)]:
                train_ppis_overlap.add(train_ppi)
    print "removing train_ppis overlapped with test clusters: %s" % len(train_ppis_overlap)
    train_ppis_clean = train_ppis_clean - train_ppis_overlap


    #kdrew: generate negative ppi set for test
    #kdrew: generate set of all proteins in training and test set
    all_test_proteins = set()
    for ppi in full_test_ppis:
        all_test_proteins = all_test_proteins.union(ppi)
    #kdrew: negative ppis are all pairs of proteins without those in the positive set
    test_neg_ppis = set()
    for cpair in it.combinations(all_test_proteins,2):
        pair = frozenset(cpair)
        if pair not in full_test_ppis:
            test_neg_ppis.add(pair)

    #kdrew: generate negative ppi set for train
    all_train_proteins = set()
    for ppi in full_train_ppis:
        all_train_proteins = all_train_proteins.union(ppi)
    #kdrew: negative ppis are all pairs of proteins without those in the positive set
    train_neg_ppis = set()
    for cpair in it.combinations(all_train_proteins,2):
        pair = frozenset(cpair)
        if pair not in full_train_ppis:
            train_neg_ppis.add(pair)


    #kdrew: remove intersection from neg ppis sets
    neg_test_ppis_clean, neg_train_ppis_clean = remove_intersection_random(test_neg_ppis, train_neg_ppis)

    #kdrew: remove intersection between neg_test_ppis_clean and full_train_ppis from neg_test_ppis_clean
    neg_test_ppis_clean = remove_intersection(neg_test_ppis_clean, full_train_ppis)

    #kdrew: remove intersection between neg_train_ppis_clean and full_test_ppis from neg_train_ppis_clean
    neg_train_ppis_clean = remove_intersection(neg_train_ppis_clean, full_test_ppis)

    #kdrew: if input_complexes is set then remove all pairs from the negative ppis
    if input_complexes != None:
        print "removing intersection with pre-merged input_complexes"
        full_input_pairs = set([frozenset(z) for y in input_complexes for z in it.combinations(y,2)])
        neg_test_ppis_clean = remove_intersection(neg_test_ppis_clean, full_input_pairs)
        neg_train_ppis_clean = remove_intersection(neg_train_ppis_clean, full_input_pairs)

    if remove_large_complexes:
        test_list_clean = [x for x in test_list_clean if len(x) <= size_threshold]
        train_list_clean = [x for x in train_list_clean if len(x) <= size_threshold]

    return test_list_clean, train_list_clean, test_ppis_clean, train_ppis_clean, neg_test_ppis_clean, neg_train_ppis_clean

#kdrew: remove overlapping complexes randomly from test and training sets
def remove_overlapping_complexes(complex_list1, complex_list2):
    #kdrew: make copies of input complex lists so they can be manipulated
    complex_list1_processed = [set(c1) for c1 in complex_list1]
    complex_list2_processed = [set(c2) for c2 in complex_list2]

    #kdrew: keep track of which list each complex came from and store all in a single list
    complex_list = [(1,c1) for c1 in complex_list1_processed]
    complex_list = complex_list + [(2,c2) for c2 in complex_list2_processed]

    #kdrew: shuffle full list of complexes
    np.random.shuffle(complex_list)

    #kdrew: for each complex in shuffled list
    for i, c1 in complex_list:
        #kdrew: if the complex is in complex list 1 and overlaps with any complexes in list 2, remove from list 1
        if i == 1 and np.any([len(c2.intersection(c1))>1 for c2 in complex_list2_processed]):
            #kdrew: remove complex from complex_list1_processed
            complex_list1_processed.remove(c1)

        #kdrew: if the complex is in complex list 2 and overlaps with any complexes in list 1, remove from list 2
        elif i == 2 and np.any([len(c2.intersection(c1))>1 for c2 in complex_list1_processed]):
            #kdrew: remove complex from complex_list2_processed
            complex_list2_processed.remove(c1)

    return complex_list1_processed, complex_list2_processed


#kdrew: removes intersection from ppi_set1
def remove_intersection(ppi_set1, ppi_set2):
    #kdrew: find intersection between ppi_set1 and ppi_set2
    intersection = list(ppi_set1.intersection(ppi_set2))

    print "removing %s" % (len(set(intersection)))
    #kdrew: remove intersection from ppis_set1
    ppi_set1_clean = ppi_set1 - set(intersection)

    return ppi_set1_clean

#kdrew: will randomly split intersection and remove half from each
def remove_intersection_random(ppi_set1, ppi_set2):
    #kdrew: find intersection between ppi_set1 and ppi_set2
    intersection = list(ppi_set1.intersection(ppi_set2))

    #kdrew: shuffle intersection
    np.random.shuffle(intersection)
    
    #kdrew: split shuffled intersection in half
    split_id = len(intersection)/2
    intersection1 = intersection[:split_id]
    intersection2 = intersection[split_id:]
    #print "intersection1"
    #print intersection1
    #print "intersection2"
    #print intersection2

    #kdrew: remove half of shuffled intersection from ppi_set1
    ppi_set1_clean = ppi_set1 - set(intersection1)
    #kdrew: remove other half of shuffled intersection from ppi_set2
    ppi_set2_clean = ppi_set2 - set(intersection2)

    return ppi_set1_clean, ppi_set2_clean



if __name__ == "__main__":
    main()


