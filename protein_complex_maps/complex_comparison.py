    
import numpy as np
import numpy.random as rand
import pandas as pd
import pickle as p
import argparse
import itertools as it
import random
import bisect
import scipy.misc as misc
from scipy.stats import hmean


import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


class ComplexComparison(object):

    def __init__(self, gold_standard=[], clusters=[], exclusion_complexes=[], samples=10000, pseudocount=1, exact=False, max_clique=None):
        #kdrew: gold_standard and clusters are a list of complexes 
        #kdrew: each complex contains a set of ids  (if passed in a list, will be converted to set)
        self.gold_standard = [set(x) for x in gold_standard]

        self.gold_standard_proteins = set()
        for x in gold_standard:
            self.gold_standard_proteins = self.gold_standard_proteins.union(x)

        self.clusters = [set(x) for x in clusters]

        self.exclusion_complexes = [set(x) for x in exclusion_complexes]

        #kdrew: dataframe where columns are clusters and rows are complexes(gold standard)
        self.intersection_table = None
        self.na_table = None

        self.clique_comparison_metric_results = None
        #kdrew: number of samples for clique comparison metric
        self.samples=samples
        self.pseudocount = pseudocount
        self.exact = exact
        self.max_clique = max_clique

    def get_gold_standard(self,):
        return self.gold_standard

    def get_exclusion_complexes(self,):
        return self.exclusion_complexes

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

    #kdrew: harmonic mean of f1-scores across clique sizes, single score for comparison of cluster predictions to gold standard complexes
    def clique_comparison_metric_grandf1score(self,mean_func = hmean):
        result_dict = self.clique_comparison_metric()
        f1score_list = [result_dict[x]['f1score'] for x in result_dict.keys()]

        grand_f1score = 0.0
        if 0.0 not in f1score_list:
            grand_f1score = mean_func(f1score_list)

        return grand_f1score


    def clique_comparison_metric_mean(self,):
        result_dict = self.clique_comparison_metric()
        precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
        recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]

        #print precision_list
        #print recall_list
        
        #print "mean precision %s mean recall %s" % (np.mean(precision_list), np.mean(recall_list))

        return {'precision_mean':np.mean(precision_list),'recall_mean':np.mean(recall_list)}

    #kdrew: calculate precision, recall and f1score for all clique sizes up to largest cluster
    def clique_comparison_metric(self, force=False):

        #kdrew: if already computed, just return the saved result
        if self.clique_comparison_metric_results != None and not force:
            return self.clique_comparison_metric_results

        #kdrew: only evaluate on proteins that are in gold standard set
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters()]
        clust_max_len = np.max(map(len,clusters))
        gs_max_len = np.max(map(len,self.get_gold_standard()))

        #kdrew: find the max cluster sizes of all clusters and gold standard complexes, then find the min between the maxes 
        #kdrew: I don't know if I like this because it will not report anything for predictions that are considered false positives but are to large to be considered
        #min_of_max_len = np.min([max_len, np.max(map(len,self.get_gold_standard()))])


        return_dict = dict()
        #print "max_len %s" % max_len

        if self.max_clique == None:
            max_len = clust_max_len
        else:
            #kdrew: if max_clique is set but the largest cluster is smaller, use the smaller value
            max_len = min(self.max_clique, clust_max_len)

        cumulative_gs_tp = 0.0
        cumulative_tp = 0.0
        cumulative_fn = 0.0
        cumulative_fp = 0.0
        for size in range(2, max_len+1):

            zero_w_pseudocount = 1.0*self.pseudocount / (self.samples + 2*self.pseudocount)

            #kdrew: shortcut, if the last two precision and recall estimates are 0.0 then there is exceedingly small chance are a larger clique size will have a precision/recall > 0.0
            #kdrew: do not calculate, just mark the remainder clique sizes as zero
            if size > 4 and return_dict[size-1]['precision'] == zero_w_pseudocount and return_dict[size-1]['recall'] == zero_w_pseudocount and return_dict[size-2]['precision'] == zero_w_pseudocount and return_dict[size-2]['recall'] == zero_w_pseudocount:
                return_dict[size] = {'precision':zero_w_pseudocount,'f1score':zero_w_pseudocount, 'recall':zero_w_pseudocount, 'cumulative_recall':zero_w_pseudocount, 'cumulative_precision':zero_w_pseudocount, 'numOfClusters':0}
                continue

            #kdrew: if clique size is greater than the max size of the gold standard, return precision 0.0 and recall 0.0 
            #kdrew: (technically recall should be NA because there are no complexes left in the gold standard but I think it makes sense to reduce it to 0.0 here)
            #kdrew: see above comment on another way of dealing with this (i.e. min_of_max_len)
            if gs_max_len < size:
                return_dict[size] = {'precision':zero_w_pseudocount,'f1score':zero_w_pseudocount, 'recall':zero_w_pseudocount, 'cumulative_recall':zero_w_pseudocount, 'cumulative_precision':zero_w_pseudocount, 'numOfClusters':0}
                continue

            if self.exact:
                result_dict = self.clique_comparison_exact(size)
            else:
                result_dict = self.clique_comparison(size)

            cumulative_gs_tp += result_dict['gs_tp']
            cumulative_tp += result_dict['tp']
            cumulative_fp += result_dict['fp']
            cumulative_fn += result_dict['fn']

            cumulative_recall = 1.0*cumulative_gs_tp / (cumulative_gs_tp + cumulative_fn)
            cumulative_precision = 1.0*cumulative_tp / (cumulative_tp + cumulative_fp)

            #print result_dict
            recall = 1.0*result_dict['gs_tp'] / (result_dict['gs_tp'] + result_dict['fn'])
            precision = 1.0*result_dict['tp'] / (result_dict['tp'] + result_dict['fp'])
            #print "clique_size: %s precision: %s recall: %s" % (size, precision, recall)

            #kdrew: calculate f1score (hmean of precision and recall) for each clique size, hmean is undefined when any value is <= 0, force it to be 0.0
            f1score = 0.0
            if precision != 0.0 and recall != 0.0:
                f1score = hmean([recall,precision])

            return_dict[size] = {'precision':precision,'recall':recall,'f1score':f1score, 'cumulative_recall':cumulative_recall, 'cumulative_precision':cumulative_precision, 'numOfClusters':result_dict['numOfClusters']}
            #print return_dict[size]

        #kdrew: save results for future retrieval
        self.clique_comparison_metric_results = return_dict
        return return_dict

    def clique_comparison_exact(self,clique_size):
        true_positives = self.pseudocount
        gs_true_positives = self.pseudocount
        false_positives = self.pseudocount
        false_negatives = self.pseudocount

        if len(self.get_exclusion_complexes()) > 0:
            raise Exclusion_Complexes_Exception("ERROR: EXCLUSION COMPLEX functionality is not implemented for exact method (def clique_comparison_exact), use estimated (def clique_comparison)")


        #kdrew: only get clusters that are larger than or equal to the clique size
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters() if len(clust & self.get_gold_standard_proteins()) >= clique_size]

        for clust in clusters:
            is_positive_list = [ np.max(map(set(group).issubset,self.get_gold_standard())) for group in it.combinations(clust, clique_size)]
            tp_curr = sum(is_positive_list)
            true_positives += tp_curr

            false_positives += len(is_positive_list) - tp_curr 

        for gs_clust in self.get_gold_standard():
            is_positive_list = [ np.max(map(set(gs_group).issubset,clusters)) for gs_group in it.combinations(gs_clust, clique_size) ]
            gs_true_positives += sum(is_positive_list)
            false_negatives += sum(np.logical_not(is_positive_list))

        return_dict = dict()
        return_dict['tp'] = true_positives
        return_dict['gs_tp'] = gs_true_positives
        return_dict['fp'] = false_positives
        return_dict['fn'] = false_negatives
        return_dict['numOfClusters'] = len(clusters)
        return return_dict



    #kdrew: calculate confusion matrix between predicted clusters and gold standard complexes for specific clique sizes
    def clique_comparison(self, clique_size):
        true_positives = self.pseudocount
        gs_true_positives = self.pseudocount
        false_positives = self.pseudocount
        false_negatives = self.pseudocount


        #kdrew: only get clusters that are larger than or equal to the clique size
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters() if len(clust & self.get_gold_standard_proteins()) >= clique_size]

        #print "clique_size: %s, #ofClusters: %s" % (clique_size, len(clusters))

        #print "clusters size: %s" % (len(clusters))

        #for clust in clusters:
        #    print "clust: %s" % clust


        #kdrew: weight each cluster by the number of combinations of size of clique overlapping with the gold standard
        #kdrew: this is so larger clusters with more possible combinations are sampled more
        weights_array = [misc.comb(len(clust & self.get_gold_standard_proteins()), clique_size) for clust in clusters ] 
        #print clique_size
        #print weights_array
        wrg = WeightedRandomGenerator( weights_array )

        #kdrew: generate set of random cliques
        random_cliques = set()
        for s in xrange(self.samples):

            #kdrew: get a random cluster weighted by the number of possible clique combinations
            clust = clusters[wrg()]
            #print clust

            #kdrew: only look at cluster ids that are in the gold standard
            clust_intersection = clust & self.get_gold_standard_proteins()
            #print clust_intersection

            #kdrew: if all proteins are outside of the gold standard, move on (unlikely to ever get selected due to weighted sampling)
            if len(clust_intersection) <= 0:
                continue

            #shuffled_l = rand.permutation(list(clust))
            shuffled_l = rand.permutation(list(clust_intersection))

            random_cliques.add(frozenset(shuffled_l[:clique_size]))

        for random_clique in random_cliques:
            #kdrew: the max just finds if any are true, could probably use .any here
            if np.max(map(random_clique.issubset,self.get_gold_standard())):
                true_positives += 1 
            else:
                if len(self.get_exclusion_complexes()) > 0 and np.max(map(random_clique.issubset,self.get_exclusion_complexes())):
                    #print "excluded clique"
                    continue
                else:
                    #print "false_positive: %s" % (' '.join(random_clique),)
                    false_positives += 1


        #kdrew: only get gold standard complexes that are larger than or equal to the clique size
        gs_clusters = [gs_clust for gs_clust in self.get_gold_standard() if len(gs_clust) >= clique_size]
        #print clique_size
        #print "gs_clusters size: %s" % (len(gs_clusters))

        #for gs_clust in gs_clusters:
            #print "gs_clust: %s" % (len(gs_clust),)
            #print "gs_clust: %s" % (gs_clust,)

        #kdrew: weight each complex by number of possible combinations of clique size for each complex
        #kdrew: NOTE: all individual cliques should be treated equally but this might be overweighting cliques that are in multiple complexes, not sure how to fix yet 
        #kdrew: this should be fixed now with the random_gs_cliques set, only evaluates a clique once
        gs_wrg = WeightedRandomGenerator( [misc.comb(len(gs_clust), clique_size) for gs_clust in gs_clusters ] )

        #kdrew: generate set of random cliques
        random_gs_cliques = set()
        for s in xrange(self.samples):

            #kdrew: get a random cluster
            gs_clust = gs_clusters[gs_wrg()]

            shuffled_l = rand.permutation(list(gs_clust))
            #print "sampled: %s" % (shuffled_l[:clique_size],)

            random_gs_cliques.add(frozenset(shuffled_l[:clique_size]))

        for random_clique in random_gs_cliques:

            #if np.max(map( set( shuffled_l[:clique_size] ).issubset, clusters )):
            if np.max(map( random_clique.issubset, clusters )):
                gs_true_positives += 1 
            else:
                false_negatives += 1

        #print "truepos: %s gs_truepos: %s falsepos: %s falseneg: %s" % (true_positives, gs_true_positives, false_positives, false_negatives)

        #assert true_positives == gs_true_positives

        return_dict = dict()
        return_dict['tp'] = true_positives
        return_dict['gs_tp'] = gs_true_positives
        return_dict['fp'] = false_positives
        return_dict['fn'] = false_negatives
        return_dict['numOfClusters'] = len(clusters)
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

class Exclusion_Complexes_Exception(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

if __name__ == "__main__":
        main()

