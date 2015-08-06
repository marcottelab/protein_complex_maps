    
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
 
import pandas as pd
import pickle as p
import argparse

import protein_complex_maps.complex_comparison as cc


def main():

    parser = argparse.ArgumentParser(description="Compare cluster predictions to gold standard complexes")
    parser.add_argument("--cluster_predictions", action="store", dest="cluster_filename", required=True, 
                                            help="Filename of cluster predictions, format one cluster per line, ids space separated")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True, 
                                            help="Filename of gold standard complexes, format one complex per line, ids space separated")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False,  default=None,
                                            help="Filename of plot")

    args = parser.parse_args()


    #fmap = open("/home/kdrew/data/bioplex/mmc4-2_tableS3_clusterMembers.geneids2uniprot","rb")
    #umap = dict()
    #for line in fmap.readlines():
    #    umap[line.split()[0]] = line.split()[1]
    #    umap[line.split()[1]] = line.split()[0]
    #
    #nhoods = p.load(open("/home/kdrew/data/bioplex/mmc4-2_tableS3_neighborhoods.p",'rb'))
    #
    #def map_id(id1):
    #    try: 
    #        return umap[id1]
    #    except KeyError:
    #        return ""
    #
    #clpred_f = open(args.cluster_filename,"wb")
    #uniprot_nhoods = []
    #for key in sorted([int(k) for k in nhoods.keys()]):
    #    clpred_f.write(' '.join(map(map_id, nhoods[str(key)])))
    #    clpred_f.write('\n')
    #    uniprot_nhoods.append( map(map_id, nhoods[str(key)]))
    #
    #clpred_f.close()

    clpred_f = open(args.cluster_filename,"rb")
    predicted_clusters = []
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())

    fcorum = open(args.gold_standard_filename,"rb")

    gold_standard_complexes = []
    for line in fcorum.readlines():
        gold_standard_complexes.append(line.split())

    cplx_compare = cc.ComplexComparison(gold_standard_complexes, predicted_clusters)
    #print cplx_compare.sensitivity(topN=10)
    #print cplx_compare.ppv(topN=10)
    #print cplx_compare.acc(topN=10)
    #print cplx_compare.sensitivity()
    #print cplx_compare.ppv()
    #print cplx_compare.acc()

    sensitivity_results = []
    ppv_results = []
    acc_results = []
    for i in range(len(predicted_clusters)):
        sensitivity_results.append(cplx_compare.sensitivity(topN=i))
        ppv_results.append(cplx_compare.ppv(topN=i))
        acc_results.append(cplx_compare.acc(topN=i))

    plt.plot(range(len(predicted_clusters)), sensitivity_results,'blue',label='Sensitivity')
    plt.plot(range(len(predicted_clusters)), ppv_results, 'red', label='PPV')
    plt.plot(range(len(predicted_clusters)), acc_results, 'grey', label='Accuracy')
    plt.legend(loc='lower right')
    plt.xlabel('Rank')

    if args.plot_filename is None:
        print "savefilename is None"
        plt.show()
    else:
        plt.savefig(args.plot_filename)


if __name__ == "__main__":
        main()



