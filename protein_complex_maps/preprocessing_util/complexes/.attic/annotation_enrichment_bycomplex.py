    
import numpy as np
import pandas as pd
import pickle as p
import argparse
import itertools as it
import csv

from scipy.misc import comb
from scipy.stats import hypergeom

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

import protein_complex_maps.features.shared_bait_feature as sbf
import protein_complex_maps.protein_util as pu

import protein_complex_maps.complexes.annotation_enrichment as ae

#kdrew: stores results of hypergeometric test
#kdrew: storing pval_corr here might be a design flaw because it decouples it from the other results which the correction is based on
class Result(object):
    def __init__(self, k=None, m=None, n=None, N=None, pval=None, pval_corr=None):
        self.k = k
        self.m = m
        self.n = n
        self.N = N
        self.pval = pval
        self.pval_corr = pval_corr

def main():

    parser = argparse.ArgumentParser(description="Enrichment analysis for each complex")
    parser.add_argument("--complexes", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes, format one cluster per line, ids space separated")
    parser.add_argument("--annotations", action="store", dest="annotation_filename", required=True, 
                                            help="Filename of annotations")
    parser.add_argument("--format", action="store", dest="file_format", required=True, 
                                            help="File format of annotations file: gsea, pairs or MIM_bliebeskind")

    parser.add_argument("--from_id", action="store", dest="from_id", required=False, default=None,
                                            help="annotation protein id type, (ENSEMBL_ID,etc), default=None (same as complexes protein id type)")
    parser.add_argument("--to_id", action="store", dest="to_id", required=False, default="ACC", 
                                            help="convert annotation ids to this type, (P_ENTREZGENEID,ACC,etc), default=ACC")
    parser.add_argument("--reviewed", action="store_true", dest="reviewed", required=False, default=False,
                                            help="map only to reviewed ids, default=False")
    parser.add_argument("--pval_correct", action="store", dest="pval_correct", required=False, default="BH",
                                            help="Type of pvalue correction, types are from R's p_adjust routine, default=BH")

    args = parser.parse_args()


    complexes = []
    complex_proteins = set()
    complex_file = open(args.complex_filename,"rb")
    for line in complex_file.readlines():
        complexes.append(line.split())
        complex_proteins = complex_proteins.union(set(line.split()))

    complex_file.close()

    annotations = ae.Annotations(filename=args.annotation_filename, file_format=args.file_format, from_id=args.from_id, to_id=args.to_id, reviewed=args.reviewed )

    stats = importr('stats')

    results = calc_enrichment(complexes, annotations)
    for i, c in enumerate(complexes):

        print "%s %s" % ( i, ' '.join(c) )

        if args.pval_correct != None:
            pvals = [results[i][term].pval for term in results[i]]
            p_adjust = stats.p_adjust(FloatVector(pvals), method = args.pval_correct)
            for j, term in enumerate(results[i]):
                results[i][term].pval_corr = p_adjust[j]

        for term in results[i]:
            print term
            print "pval_corr: %s, pval: %s, k: %s, m: %s, n: %s, N: %s" % (results[i][term].pval_corr, results[i][term].pval, results[i][term].k,  results[i][term].m,  results[i][term].n,  results[i][term].N)  
            print " "


#kdrew: calculate enrichment for protein pairs with shared annotations
#kdrew: protein ids in both complexes and protein_annotaiton_dict must be of the same type
def calc_enrichment(complexes, annotations):

    protein_annotation_dict = annotations.get_protein_annotation_dict()
    annotation_dict = annotations.get_annotation_dict()

    #kdrew: generate set of proteins in complexes
    complex_proteins = set([p for c in complexes for p in c])
    #kdrew: generate set of annotated proteins
    annotated_proteins = set(protein_annotation_dict.keys())

    #print complex_proteins
    #print annotated_proteins
    annotated_and_complex_proteins = complex_proteins.intersection(annotated_proteins)


    #kdrew: calculate N, m, n and k
    #kdrew: where N is the # of total complex proteins with annotations
    print "number of annotated proteins: %s" % len(annotated_and_complex_proteins)
    N = len(annotated_and_complex_proteins)

    complex_results = dict()
    for i, complex1 in enumerate(complexes):
        #print complex1
        term_results = dict()
        for term in annotation_dict:

            #kdrew: m is # of proteins annotated with term
            m = len(set(annotation_dict[term]).intersection(annotated_and_complex_proteins))

            #kdrew: n is # of annotated complex proteins
            n = len(set(complex1).intersection(annotated_and_complex_proteins))

            #kdrew: k is # of complex proteins with term
            k = len(set(complex1).intersection(set(annotation_dict[term])))

            #pval = sbf.pval(k,n,m,N, adhoc=True)
            #pval = sbf.pval(k,n,m,N, adhoc=False)
            pval = hypergeom(N,m,n).pmf(k)

            #kdrew: only store the pvalues from terms that have proteins with its annotation and are in the complex with the term
            #kdrew: gprofiler has a description of how to choose the number of tests for correction: http://biit.cs.ut.ee/gprofiler/help.cgi?help_id=7
            if m > 0 and k > 0:
                #print term
                #print "k: %s, n: %s, m: %s, N: %s" % (k,n,m,N)
                #print "pval: %s" % pval
                term_results[term] = Result(k=k,n=n,m=m,N=N,pval=pval)

        complex_results[i] = term_results

    return  complex_results

if __name__ == "__main__":
        main()


