
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it
import operator

def main():

    parser = argparse.ArgumentParser(description="Split gold standard complexes into test and training")
    parser.add_argument("--input_complexes", action="store", dest="complex_filename", required=True, 
                                    help="Filename of gold standard complexes")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", required=False, default=None,
                                    help="Random seed for shuffling, default=None")
    args = parser.parse_args()

    if args.random_seed != None:
        print "setting random seed"
        np.random.seed(args.random_seed)

    #kdrew: origin of this code is from ipython notebook: corum_test_train_revisit, heavily refactored and modified

    complexes = []
    f = open(args.complex_filename,"rb")
    for line in f.readlines():
        complexes.append(line.split())

    test_list_clean, train_list_clean, test_ppis_clean, train_ppis_clean, neg_test_ppis_clean, neg_train_ppis_clean = split_complexes(complexes)

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
def split_complexes(complexes): 

    #kdrew: randomly shuffle complexes 
    np.random.shuffle(complexes)
    split_id = len(complexes)/2
    test_list = complexes[:split_id]
    train_list = complexes[split_id:]

    test_list_clean, train_list_clean = remove_overlapping_complexes(test_list, train_list)

    #print "test_list"
    #print test_list
    #print "train_list"
    #print train_list
    
    #kdrew: for complexes in test_list, generate all pairs
    test_ppis = set()
    for clust in test_list:
        for pair in it.combinations( set( clust ), 2):
            test_ppis.add(frozenset(pair))

    #print "test_ppis: "
    #print test_ppis

    #kdrew: for complexes in train_list, generate all pairs
    train_ppis = set()
    for clust in train_list:
        for pair in it.combinations( set( clust ), 2):
            train_ppis.add(frozenset(pair))

    #print "train_ppis: "
    #print train_ppis

    #kdrew: remove intersection from test and train ppi sets
    test_ppis_clean, train_ppis_clean = remove_intersection_random(test_ppis, train_ppis)

    #print "test_ppis_clean: "
    #print test_ppis_clean

    #print "train_ppis_clean: "
    #print train_ppis_clean
    
    

    #kdrew: generate negative ppi set for test
    #kdrew: generate set of all proteins in training and test set
    all_test_proteins = set()
    for ppi in test_ppis:
        all_test_proteins = all_test_proteins.union(ppi)
    #kdrew: negative ppis are all pairs of proteins without those in the positive set
    test_neg_ppis = set()
    for cpair in it.combinations(all_test_proteins,2):
        pair = frozenset(cpair)
        if pair not in test_ppis:
            test_neg_ppis.add(pair)

    #kdrew: generate negative ppi set for train
    all_train_proteins = set()
    for ppi in train_ppis:
        all_train_proteins = all_train_proteins.union(ppi)
    #kdrew: negative ppis are all pairs of proteins without those in the positive set
    train_neg_ppis = set()
    for cpair in it.combinations(all_train_proteins,2):
        pair = frozenset(cpair)
        if pair not in train_ppis:
            train_neg_ppis.add(pair)


    #kdrew: remove intersection from neg ppis sets
    neg_test_ppis_clean, neg_train_ppis_clean = remove_intersection_random(test_neg_ppis, train_neg_ppis)

    #kdrew: remove intersection between neg_test_ppis_clean and train_ppis from neg_test_ppis_clean
    neg_test_ppis_clean = remove_intersection(neg_test_ppis_clean, train_ppis)

    #kdrew: remove intersection between neg_train_ppis_clean and test_ppis from neg_train_ppis_clean
    neg_train_ppis_clean = remove_intersection(neg_train_ppis_clean, test_ppis)


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


