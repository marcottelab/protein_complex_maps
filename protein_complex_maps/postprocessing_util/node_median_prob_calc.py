
import numpy as np
import argparse
import itertools as it

import pandas as pd

def main():

    parser = argparse.ArgumentParser(description="Calculate median probability of each protein for each complex")
    parser.add_argument("--complexes", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes, format one cluster per line, ids space separated")
    parser.add_argument("--pairwise_filename", action="store", dest="pairwise_filename", required=True,
                                            help="Filename of pairwise file, (ie. one line per ppi, third column contains score)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=False, default=None, 
                                            help="Filename to output median probability scores")
    args = parser.parse_args()


    complexes = []
    f = open(args.complex_filename, "r")
    for clustid, line in enumerate(f.readlines()):
        #c = ["%s_%s" % (clustid, x) for x in line.split()]
        complexes.append(set(line.split()))

    pairwise_interactions = dict()
    f2 = open(args.pairwise_filename, "r")
    for line in f2.readlines():
        try:
            pairwise_interactions[frozenset([line.split()[0],line.split()[1]])] = line.split()[2]
        except IndexError:
            continue



    prot_zscores = dict()
    prot_medians = dict()
    for clust_id, c in enumerate(complexes):
        #print("clust_id: %s" % clust_id)
        clust_prot_dict = dict()
        clust_list = []
        for prot_pair in it.combinations(c,2):
            try:
                #print("prot_pair: %s = %s" % (prot_pair, pairwise_interactions[frozenset(prot_pair)]))
                score = float(pairwise_interactions[frozenset(prot_pair)])
                clust_list.append(score)
                #kdrew: add score of first protein
                try:
                    clust_prot_dict[prot_pair[0]].append(score)
                except KeyError:
                    clust_prot_dict[prot_pair[0]] = [score]
                try:
                    clust_prot_dict[prot_pair[1]].append(score)
                except KeyError:
                    clust_prot_dict[prot_pair[1]] = [score]
            except KeyError:
                #print("prot_pair: %s = %s" % (prot_pair, 0.0))
                clust_list.append(0.0)
                try:
                    clust_prot_dict[prot_pair[0]].append(0.0)
                except KeyError:
                    clust_prot_dict[prot_pair[0]] = [0.0]
                try:
                    clust_prot_dict[prot_pair[1]].append(0.0)
                except KeyError:
                    clust_prot_dict[prot_pair[1]] = [0.0]

        #print(clust_prot_dict)
        #print(clust_list)
        clust_medians = dict()
        for prot in clust_prot_dict:
            prot_median = np.median(clust_prot_dict[prot])
            #print("%s : %s" % (prot, prot_median ))
            clust_medians["%s_%s" % (clust_id,prot)] = prot_median

        #print(clust_medians)
        prot_medians.update(clust_medians)
        #print(list(clust_medians.values()))
        clust_medians_mean = np.mean(list(clust_medians.values()))
        clust_medians_std = np.std(list(clust_medians.values()))
        for prot in clust_medians:
            zscore = (clust_medians[prot] - clust_medians_mean)/clust_medians_std
            prot_zscores[prot] = zscore

    #print(prot_zscores)
    #print(prot_medians)
    m_df = pd.DataFrame().from_dict(prot_medians, orient='index',columns=['median_complex_edge_score'])
    z_df = pd.DataFrame().from_dict(prot_zscores, orient='index',columns=['zscore_complex_edge_score'])

    m_df.join(z_df).to_csv(args.output_filename)


if __name__ == "__main__":
        main()

