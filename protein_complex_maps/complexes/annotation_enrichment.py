    
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn3_circles, venn2
import pandas as pd
import pickle as p
import argparse
import itertools as it
import csv

from scipy.misc import comb
from scipy.stats import hypergeom

import protein_complex_maps.features.shared_bait_feature as sbf
import protein_complex_maps.protein_util as pu




class Annotations(object):

    def __init__(self, filename=None, file_format=None, protein_annotation_dict=None, annotation_dict=None, from_id=None, to_id=None, reviewed=False):
        self.__filename = filename
        self.__file_format = file_format
        self.__protein_annotation_dict = protein_annotation_dict
        self.__annotation_dict = annotation_dict

        if self.__filename != None:
            if self.__file_format == "gsea":
                self.__annotation_dict = dict()
                annotation_file = open(self.__filename, "rb")
                for line in annotation_file.readlines():
                    term = line.split()[0]
                    term_url = line.split()[1]
                    protein_ids = line.split()[2:]
                    self.__annotation_dict[term] = protein_ids

            elif self.__file_format == "pairs":
                self.__protein_annotation_dict = dict()
                annotation_file = open(self.__filename, "rb")
                for line in annotation_file.readlines():
                    prot_id = line.split()[0]
                    term = line.split()[1]

                    try:
                        protein_annotation_dict[prot_id].append(term)
                    except KeyError:
                        protein_annotation_dict[prot_id] = [term]
            elif self.__file_format == "MIM_bliebeskind":
                self.__annotation_dict = dict()
                annotation_file = open(self.__filename, "rb")
                reader = csv.reader(annotation_file)
                for row in reader:
                    if "group" in row[0] and "percentLECA" in row[1]:
                        continue
                    else:
                        protein_ids = row[5].split('|')
                        term = row[7].strip()
                        #print term
                        #print protein_ids
                        self.__annotation_dict[term] = protein_ids

            else:
                print "no file_format"

        if from_id != None and to_id != None:
            if self.__protein_annotation_dict != None:

                inputID2ACC_map = pu.map_protein_ids(self.__protein_annotation_dict.keys(), from_id, "ACC", reviewed=reviewed)
                flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
                ACC2outputID_map = pu.map_protein_ids(flatten_list, "ACC", to_id, reviewed=reviewed)

                mapped_protein_annotation_dict = dict()
                for prot_id in self.__protein_annotation_dict:
                    mapped_id = ACC2outputID_map[inputID2ACC_map[prot_id]]
                    mapped_protein_annotation_dict[mapped_id] = self.__protein_annotation_dict[prot_id]

                self.__protein_annotation_dict = mapped_protein_annotation_dict

            elif self.__annotation_dict != None:
                print self.__annotation_dict.values()
                protein_ids = list(set([x for l in self.__annotation_dict.values() for x in l]))
                inputID2ACC_map = pu.map_protein_ids(protein_ids, from_id, "ACC", reviewed=reviewed)
                flatten_list = [item for sublist in inputID2ACC_map.values() for item in sublist]
                ACC2outputID_map = pu.map_protein_ids(flatten_list, "ACC", to_id, reviewed=reviewed)

                mapped_annotation_dict = dict()
                for term in self.__annotation_dict:
                    for prot_id in self.__annotation_dict[term]:
                        for acc in inputID2ACC_map[prot_id]:
                            if len(ACC2outputID_map[acc]) == 0:
                                continue
                            else:
                                mapped_id = ACC2outputID_map[acc][0]
                                break
                        try:
                            mapped_annotation_dict[term].append(mapped_id)
                        except KeyError:
                            mapped_annotation_dict[term] = [mapped_id]

                self.__annotation_dict = mapped_annotation_dict


        #kdrew: populate protein_annotation_dict
        if self.__annotation_dict != None and self.__protein_annotation_dict == None:
            self.__protein_annotation_dict = dict()
            for term in self.__annotation_dict:
                for protid in self.__annotation_dict[term]:
                    try: 
                        self.__protein_annotation_dict[protid].append(term)
                    except KeyError:
                        self.__protein_annotation_dict[protid] = [term]

        #kdrew: populate annotation_dict
        elif self.__annotation_dict == None and self.__protein_annotation_dict != None:
            self.__annotation_dict = dict()
            for protid in self.__protein_annotation_dict:
                for term in self.__protein_annotation_dict[protid]:
                    try:
                        self.__annotation_dict[term].append(protid)
                    except KeyError:
                        self.__annotation_dict[term] = [protid]


    def get_annotation_dict(self,):
        return self.__annotation_dict

    def get_protein_annotation_dict(self,):
        return self.__protein_annotation_dict


def main():

    parser = argparse.ArgumentParser(description="Enrichment analysis of annotated protein pairs")
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
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                        help="Filename of output plot")

    args = parser.parse_args()


    complexes = []
    complex_proteins = set()
    complex_file = open(args.complex_filename,"rb")
    for line in complex_file.readlines():
        complexes.append(line.split())
        complex_proteins = complex_proteins.union(set(line.split()))

    complex_file.close()

    annotations = Annotations(filename=args.annotation_filename, file_format=args.file_format, from_id=args.from_id, to_id=args.to_id, reviewed=args.reviewed )

    results_dict = calc_enrichment(complexes, annotations)
    print results_dict['pval']

    if args.plot_filename != None:
        plot_venn(results_dict, args.plot_filename)

def plot_venn(results_dict, plot_filename, labels=['m','n']):
    m = results_dict['m']
    n = results_dict['n']
    k = results_dict['k']
    N = results_dict['N']
    pval = results_dict['pval']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.text(2, 6, r'an equation: $E=mc^2$', fontsize=15)

    mpl.rcParams.update({'font.size': 8})

    v = venn2(subsets = (m-k, n-k, k), set_labels=labels)
    plt.annotate('k, pval: %s' % (pval,), xy=v.get_label_by_id('11').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                         ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                                      arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    #plt.annotate('m', xy=v.get_label_by_id('10').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    #                     ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #                                  arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    #plt.annotate('n', xy=v.get_label_by_id('01').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    #                     ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #                                  arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))

    ax.text(-0.6,0.6, 'N: %s' % (int(N),)) 

    #ax.add_patch(
    #    mpatches.Rectangle(
    #        (-0.65, -0.65),   # (x,y)
    #        1.30,          # width
    #        1.30,          # height
    #        fill=False,
    #    )
    #)
    #
    #ax.autoscale()

    plt.savefig(plot_filename)

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
    #kdrew: where N is the # of possible annotated pairs
    print "number of annotated proteins: %s" % len(annotated_and_complex_proteins)
    N = comb(len(annotated_and_complex_proteins),2)

    #kdrew: m is # of possible pairs with shared annotation
    #kdrew: set of pairs of proteins in complex map that share annotations
    shared_annotation_pairs = set()
    for term in annotation_dict:
        #shared_annotation_pairs = shared_annotation_pairs.union([frozenset(pair) for pair in it.combinations(annotation_dict[term],2)])
        #kdrew: set of proteins in complexes annotated with term
        annotated_and_complex_proteins_for_term = complex_proteins.intersection(set(annotation_dict[term]))
        #kdrew: set of all possible pairs that share this term
        annotation_pairs_for_term = [frozenset(pair) for pair in it.combinations(annotated_and_complex_proteins_for_term, 2)]
        #kdrew: add unique pairs to total set
        shared_annotation_pairs = shared_annotation_pairs.union(annotation_pairs_for_term)

    m = len(shared_annotation_pairs)

    #kdrew: n is # of annotated interaction partners
    n = 0
    #kdrew: set of pairs of annotated proteins that are in complex together
    annotated_complex_pairs = set()
    for c in complexes:
        #kdrew: annotated proteins in complex
        annotated_c = annotated_and_complex_proteins.intersection(set(c))
        #kdrew: set of possible pairs in annotated complex
        annotated_complex_pairs_for_complex = [frozenset(pair) for pair in it.combinations(annotated_c, 2)]
        #kdrew: add unique pairs to total set
        annotated_complex_pairs = annotated_complex_pairs.union(annotated_complex_pairs_for_complex)

    n = len(annotated_complex_pairs)

    #kdrew: k is # of interaction partners with shared annotation
    k = len(shared_annotation_pairs.intersection(annotated_complex_pairs))
    #k = 0
    #for c in complexes:
    #    for pair in it.combinations(c,2):
    #        try:
    #            prot1_annotations = set(protein_annotation_dict[pair[0]])
    #            prot2_annotations = set(protein_annotation_dict[pair[1]])
    #            #kdrew: any overlap between protein 1 and protein 2 annotations
    #            if 0 < len(prot1_annotations.intersection(prot2_annotations)):
    #                k += 1
    #        except KeyError:
    #            continue
    

    #pval = sbf.pval(k,n,m,N, adhoc=True)
    #pval = sbf.pval(k,n,m,N, adhoc=False)
    pval = 0.0
    for i in range(k,int(min(m,n)+1)):
        pval_tmp = hypergeom(N,m,n).pmf(i)
        pval += pval_tmp
        print i, pval_tmp

    print "k: %s, n: %s, m: %s, N: %s" % (k,n,m,N)
    print "pval: %s" % pval

    return_dict = {'k':k,'n':n,'m':m,'N':N,'pval':pval}

    return return_dict


if __name__ == "__main__":
        main()


