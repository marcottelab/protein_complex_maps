    
import numpy as np
import pandas as pd
import pickle as p
import argparse

import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


class ComplexComparison(object):

    def __init__(self, gold_standard=[], clusters=[]):
        #kdrew: gold_standard and clusters are a list of complexes 
        #kdrew: each complex contains a set of ids  (if passed in a list, will be converted to set)
        self.gold_standard = [set(x) for x in gold_standard]
        self.clusters = [set(x) for x in clusters]

        #kdrew: dataframe where columns are clusters and rows are complexes(gold standard)
        self.intersection_table = None
        self.na_table = None

    def get_gold_standard(self,):
        return self.gold_standard

    def get_clusters(self,):
        return self.clusters

    def get_na_table(self,):
        if self.na_table is None:
            self.generate_na_table()
        return self.na_table

    def get_intersection_table(self,):
        if self.intersection_table is None:
            self.generate_intersection_table()
        return self.intersection_table    

    #kdrew: different metrics for comparing complexes and clusters

    #kdrew: convience wrapper function
    def mmr(self,topN=None):
        return self.max_matching_ratio(topN)

    #kdrew: method described in Yang et al. BMC Medical Genomics Integrating PPI datasets with the PPI data from biomedical literature for protein complex detection (2014)
    def max_matching_ratio(self,topN=None):
        max_df = self.max_matching_ratio_distribution(topN=topN)
        sum_max_df = max_df.sum()

        num_of_gold_complexes = len(self.get_gold_standard())

        mmr = sum_max_df / num_of_gold_complexes
        return mmr

    def max_matching_ratio_distribution(self,topN=None):
        df = self.get_na_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        #print df
        max_df = df.max(axis=1)

        return max_df


    #kdrew: method from Brohee and Helden BMC Bioinformatics, http://www.biomedcentral.com/content/pdf/1471-2105-7-488.pdf
    #kdrew: intersection/size_of_gold = percent coverage of given complex
    #kdrew: find max
    #kdrew: size_of_gold * max_percent_coverage  / size_of_gold
    def sensitivity_distribution(self, topN=None):
        df = self.get_intersection_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        sn_df = df.divide(map(len,self.get_gold_standard()), axis=0)
        #print sn_df.max()
        sn_co = sn_df.max(axis=1)
        return sn_co

    def sensitivity(self, topN=None):
        sn_co = self.sensitivity_distribution(topN=topN)
        #print sn_co.max()
        goldstd_sizes = map(len,self.get_gold_standard())
        Sn = 1.0*sum(sn_co * goldstd_sizes) / sum(goldstd_sizes)
        return Sn
        
    #kdrew: method from Brohee and Helden BMC Bioinformatics, http://www.biomedcentral.com/content/pdf/1471-2105-7-488.pdf
    def ppv_distribution(self, topN=None):
        df = self.get_intersection_table()
        #kdrew: if want to do a range of clusters (topN) we probably want to slice df here
        if topN != None:
            df = df.ix[:,:topN]
        ppv_df = df.divide(1.0*df.sum())
        ppv_df = ppv_df.fillna(0.0)
        PPV_cl = ppv_df.max()
        return PPV_cl

    def ppv(self,topN=None):
        df = self.get_intersection_table()
        PPV_cl = self.ppv_distribution(topN=topN)
        try:
            PPV = 1.0*sum(PPV_cl * df.sum()) / sum(df.sum())
        except ZeroDivisionError:
            PPV = None
        return PPV

    #kdrew: method from Brohee and Helden BMC Bioinformatics, http://www.biomedcentral.com/content/pdf/1471-2105-7-488.pdf
    def acc(self, topN=None):
        Sn = self.sensitivity(topN)
        PPV = self.ppv(topN)
        #kdrew: geometric mean
        Acc = (Sn * PPV) ** (1.0/2)
        return Acc
        
        
    def generate_intersection_table(self,):
        rows_list = []
        for cc in self.get_gold_standard():
            d = dict()
            for i, clst in enumerate(self.get_clusters()):
                d[i] = 1.0*len(set.intersection(set(clst),set(cc)))
            rows_list.append(d)
        self.intersection_table = pd.DataFrame(rows_list)

    #kdrew: method described in Yang et al. BMC Medical Genomics Integrating PPI datasets with the PPI data from biomedical literature for protein complex detection (2014)
    #kdrew: this is also the Bader Score from:  http://www.biomedcentral.com/1471-2105/4/2 
    #kdrew: NA = Neighborhood Affinity
    def generate_na_table(self,):
        rows_list = []
        for cc in self.get_gold_standard():
            d = dict()
            for i, clst in enumerate(self.get_clusters()):
                numerator = (1.0*len(set.intersection(set(clst),set(cc))))**2
                denominator = 1.0 * (len(set(clst)) * len(set(cc)))

                d[i] = numerator / denominator 
            rows_list.append(d)
        self.na_table = pd.DataFrame(rows_list)

def main():

    parser = argparse.ArgumentParser(description="Compare cluster predictions to gold standard complexes")
    parser.add_argument("--cluster_predictions", action="store", dest="cluster_filename", required=True, 
                                            help="Filename of cluster predictions, format one cluster per line, ids space separated")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True, 
                                            help="Filename of gold standard complexes, format one complex per line, ids space separated")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename for plotting histograms of metrics")

    args = parser.parse_args()

    gold_standard_complexes = []
    gold_file = open(args.gold_standard_filename,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())

    predicted_clusters = []
    clpred_f = open(args.cluster_filename,"rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())

    cplx_compare = ComplexComparison(gold_standard_complexes, predicted_clusters)
    print "Sensitivity: %s" % cplx_compare.sensitivity()
    print "PPV: %s" % cplx_compare.ppv()
    print "ACC: %s" % cplx_compare.acc()
    print "MMR: %s" % cplx_compare.mmr()


    if args.plot_filename != None:
        f, subplots = plt.subplots(3)
        subplots[0].hist(cplx_compare.max_matching_ratio_distribution())
        subplots[0].set_title('MMR')
        subplots[1].hist(cplx_compare.sensitivity_distribution())
        subplots[1].set_title('Sensitivity')
        subplots[2].hist(cplx_compare.ppv_distribution())
        subplots[2].set_title('PPV')
        plt.savefig(args.plot_filename)


if __name__ == "__main__":
        main()

