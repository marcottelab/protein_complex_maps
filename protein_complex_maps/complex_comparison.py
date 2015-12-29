    
import numpy as np
import numpy.random as rand
import pandas as pd
import pickle as p
import argparse
import itertools as it
import random
import bisect
import scipy.misc as misc

import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


class ComplexComparison(object):

    def __init__(self, gold_standard=[], clusters=[]):
        #kdrew: gold_standard and clusters are a list of complexes 
        #kdrew: each complex contains a set of ids  (if passed in a list, will be converted to set)
        self.gold_standard = [set(x) for x in gold_standard]

        self.gold_standard_proteins = set()
        for x in gold_standard:
            self.gold_standard_proteins = self.gold_standard_proteins.union(x)

        self.clusters = [set(x) for x in clusters]

        #kdrew: dataframe where columns are clusters and rows are complexes(gold standard)
        self.intersection_table = None
        self.na_table = None

    def get_gold_standard(self,):
        return self.gold_standard

    def get_gold_standard_proteins(self,):
        return self.gold_standard_proteins

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

    def clique_comparison_metric_mean(self,):
        result_dict = self.clique_comparison_metric()
        precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
        recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]

        #print precision_list
        #print recall_list
        
        #print "mean precision %s mean recall %s" % (np.mean(precision_list), np.mean(recall_list))

        return {'precision_mean':np.mean(precision_list),'recall_mean':np.mean(recall_list)}

    def clique_comparison_metric(self,):
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters()]
        max_len = np.max(map(len,clusters))
        return_dict = dict()
        #print "max_len %s" % max_len
        for size in range(2, max_len+1):

            #kdrew: shortcut, if the last two precision and recall estimates are 0.0 then there is exceedingly small chance are a larger clique size will have a precision/recall > 0.0
            #kdrew: do not calculate, just mark the remainder clique sizes as zero
            if size > 4 and return_dict[size-1]['precision'] == 0.0 and return_dict[size-1]['recall'] == 0.0 and return_dict[size-2]['precision'] == 0.0 and return_dict[size-2]['recall'] == 0.0:
                return_dict[size] = {'precision':0.0,'recall':0.0}
                continue

            result_dict = self.clique_comparison(size)
            #print result_dict
            recall = 1.0*result_dict['gs_tp'] / (result_dict['gs_tp'] + result_dict['fn'])
            precision = 1.0*result_dict['tp'] / (result_dict['tp'] + result_dict['fp'])
            #print "clique_size: %s precision: %s recall: %s" % (size, precision, recall)
            return_dict[size] = {'precision':precision,'recall':recall}

        return return_dict


    #kdrew: calculate confusion matrix between predicted clusters and gold standard complexes for specific clique sizes
    def clique_comparison(self, clique_size, samples=10000):
        true_positives = 0
        gs_true_positives = 0
        false_positives = 0
        false_negatives = 0

        #kdrew: only get clusters that are larger than or equal to the clique size
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters() if len(clust & self.get_gold_standard_proteins()) >= clique_size]

        #print "clusters size: %s" % (len(clusters))

        #for clust in clusters:
        #    print "clust: %s" % clust


        #kdrew: weight each cluster by the length of its overlap with the gold standard
        wrg = WeightedRandomGenerator( [misc.comb(len(clust & self.get_gold_standard_proteins()), clique_size) for clust in clusters ] )

        for s in xrange(samples):

            #kdrew: get a random cluster
            clust = clusters[wrg()]
            #print clust

            #kdrew: only look at cluster ids that are in the gold standard
            clust_intersection = clust & self.get_gold_standard_proteins()
            #print clust_intersection

            #kdrew: if all proteins are outside of the gold standard, move on (unlikely to ever get selected due to weighted sampling)
            if len(clust_intersection) <= 0:
                continue

            shuffled_l = rand.permutation(list(clust))

            if np.max(map(set(shuffled_l[:clique_size]).issubset,self.get_gold_standard())):
                true_positives += 1 
            else:
                false_positives +=1


        #kdrew: only get gold standard complexes that are larger than or equal to the clique size
        gs_clusters = [gs_clust for gs_clust in self.get_gold_standard() if len(gs_clust) >= clique_size]
        #print "gs_clusters size: %s" % (len(gs_clusters))

        #for gs_clust in gs_clusters:
            #print "gs_clust: %s" % (len(gs_clust),)
            #print "gs_clust: %s" % (gs_clust,)

        #kdrew: weight each complex by size of complex
        gs_wrg = WeightedRandomGenerator( [misc.comb(len(gs_clust), clique_size) for gs_clust in gs_clusters ] )

        for s in xrange(samples):

            #kdrew: get a random cluster
            gs_clust = gs_clusters[gs_wrg()]

            shuffled_l = rand.permutation(list(gs_clust))
            #print "sampled: %s" % (shuffled_l[:clique_size],)

            #if np.max(map(set(shuffled_l[:clique_size]).issubset,self.get_clusters())).any():
            if np.max(map( set( shuffled_l[:clique_size] ).issubset, clusters )):
                gs_true_positives += 1 
            else:
                false_negatives+=1

        #print "truepos: %s gs_truepos: %s falsepos: %s falseneg: %s" % (true_positives, gs_true_positives, false_positives, false_negatives)

        #assert true_positives == gs_true_positives

        return_dict = dict()
        return_dict['tp'] = true_positives
        return_dict['gs_tp'] = gs_true_positives
        return_dict['fp'] = false_positives
        return_dict['fn'] = false_negatives
        return return_dict


        
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


#kdrew: return chunks of size n from list l
def chunks(l,n):
    n = max(1, n)
    #kdrew: generator that slices list into n chunks, 
    #kdrew: if the last slice is smaller than n, ignore
    return [l[i:i + n] for i in xrange(0, len(l), n) if i+n< len(l)]

#kdrew: generator to return m number of random combinations of size n from list l
def rand_combinations(l,n,m):
    #kdrew: determine the number of shuffling needed
    try:
        i = int(np.ceil(1.0*m/(len(l)/n)))
        for ii in xrange(i):
            shuffled_l = rand.permutation(l)
            for ch in chunks(shuffled_l,n):
                yield ch
    except:
        return

#kdrew: "borrowed" code from Eli Bendersky for generating fast access to weighted lists
#http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
class WeightedRandomGenerator(object):
    def __init__(self, weights):
        self.totals = []
        running_total = 0

        for w in weights:
            running_total += w
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()



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
    ccmm = cplx_compare.clique_comparison_metric_mean()
    print "Clique Precision Mean: %s Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])


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

