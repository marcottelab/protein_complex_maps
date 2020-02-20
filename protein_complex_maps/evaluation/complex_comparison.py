    
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

    def __init__(self, gold_standard=[], clusters=[], exclusion_complexes=[], samples=10000, pseudocount=1, exact=False, max_clique=None, remove_non_gold_standard_proteins=False, normalize_by_combinations=False):
        #kdrew: gold_standard and clusters are a list of complexes 
        #kdrew: each complex contains a set of ids  (if passed in a list, will be converted to set)
        self.gold_standard = [set(x) for x in gold_standard]

        self.gold_standard_proteins = set()
        for x in gold_standard:
            self.gold_standard_proteins = self.gold_standard_proteins.union(x)

        self.cluster_proteins = set()
        self.clusters = [set(x) for x in clusters]
        for x in clusters:
            self.cluster_proteins = self.cluster_proteins.union(x)

        #kdrew: test if there are more than one protein shared between gold standard and predicted clusters
        if len(self.cluster_proteins.intersection(self.gold_standard_proteins)) <= 1:
            raise Gold_Standard_Overlap_Exception("ERROR: Gold Standard Overlap Warning: no pairs in clusters overlap with gold standard")

        if remove_non_gold_standard_proteins:
            self.remove_non_gold_standard_proteins()

        self.normalize_by_combinations = normalize_by_combinations
        

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

    #kdrew: this function removes proteins from any cluster which is not in the gold standard
    def remove_non_gold_standard_proteins(self,):
        filtered_clusters = []
        for c in self.clusters:       
            c_filt = c.intersection(self.gold_standard_proteins)
            if len(c_filt) > 1:
                filtered_clusters.append(c_filt)

        self.clusters = filtered_clusters
        return


    def get_gold_standard(self,):
        return self.gold_standard

    def get_exclusion_complexes(self,):
        return self.exclusion_complexes

    def get_gold_standard_proteins(self,):
        return self.gold_standard_proteins

    def get_clusters(self,):
        return self.clusters

    def get_cluster_proteins(self,):
        return self.cluster_proteins

    def get_na_table(self,):
        if self.na_table is None:
            self.generate_na_table()
        return self.na_table

    def get_intersection_table(self,):
        if self.intersection_table is None:
            self.generate_intersection_table()
        return self.intersection_table    

    #kdrew: different metrics for comparing complexes and clusters

    #kdrew: Song and Singh Bioinformatics 2009
    #kdrew: similar to mmr_pwmmr_hmean but weights each complex evaluation by the size of the complex
    def precision_recall_product(self,topN=None):
        pm = self.precision_measure(topN=topN)
        rm = self.recall_measure(topN=topN)
        return hmean([pm,rm])

    def precision_measure(self,topN=None):
        max_df = self.prediction_wise_max_matching_ratio_distribution(topN=topN)
        print "precision_measure: %s" % max_df
        lengths = [len(x) for x in self.get_clusters()]
        print "precision_measure: %s" % lengths
        sum_weighted_max = np.sum(lengths*max_df)
        print "precision_measure: %s" % sum_weighted_max

        precision_measure = sum_weighted_max / np.sum(lengths)
        return precision_measure

    def recall_measure(self,topN=None):
        max_df = self.max_matching_ratio_distribution(topN=topN)
        print "recall_measure: %s" % max_df
        lengths = [len(x) for x in self.get_gold_standard()]
        print "recall_measure: %s" % lengths
        sum_weighted_max = np.sum(lengths*max_df)
        print "recall_measure: %s" % sum_weighted_max

        recall_measure = sum_weighted_max / np.sum(lengths)
        return recall_measure

    def mmr_pwmmr_hmean(self,topN=None):
        mmr = self.mmr()
        pwmmr = self.pwmmr()
        return hmean([mmr,pwmmr])

    #kdrew: convience wrapper function
    def mmr(self,topN=None):
        return self.max_matching_ratio(topN)

    #kdrew: convience wrapper function
    def pwmmr(self,topN=None):
        return self.prediction_wise_max_matching_ratio(topN)

    #kdrew: method described in Yang et al. BMC Medical Genomics Integrating PPI datasets with the PPI data from biomedical literature for protein complex detection (2014)
    def max_matching_ratio(self,topN=None):
        max_df = self.max_matching_ratio_distribution(topN=topN)
        sum_max_df = max_df.sum()

        num_of_gold_complexes = len(self.get_gold_standard())

        mmr = sum_max_df / num_of_gold_complexes
        return mmr

    #kdrew: calculate mmr from the view of predicted complexes
    def prediction_wise_max_matching_ratio(self,topN=None):
        max_df = self.prediction_wise_max_matching_ratio_distribution(topN=topN)
        sum_max_df = max_df.sum()

        num_of_clusters = len(self.get_clusters())

        pwmmr = sum_max_df / num_of_clusters
        return pwmmr


    def max_matching_ratio_distribution(self,topN=None):
        df = self.get_na_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        #print df
        max_df = df.max(axis=1)

        return max_df

    #kdrew: MMR calculates the max for each standard complex, 
    #kdrew: here we calculate the max for each predicted complex
    def prediction_wise_max_matching_ratio_distribution(self,topN=None):
        df = self.get_na_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        #print df
        max_df = df.max(axis=0)

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

    def clique_comparison_metric_mean(self, weighted=False):
        result_dict = self.clique_comparison_metric()
        precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
        recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]

        #print precision_list
        #print recall_list
        
        #print "mean precision %s mean recall %s" % (np.mean(precision_list), np.mean(recall_list))

        if weighted:
            sum_complexes = sum([result_dict[x]['numOfClusters'] for x in result_dict.keys()])
            weights_list = [1.0*result_dict[x]['numOfClusters']/sum_complexes for x in result_dict.keys()]
            precision_mean = np.average(precision_list, weights=weights_list)
            recall_mean = np.average(recall_list, weights=weights_list)
        else:
            precision_mean = np.mean(precision_list)
            recall_mean = np.mean(recall_list)

        return {'precision_mean':precision_mean,'recall_mean':recall_mean}

    def clique_comparision_metric_weighted_precision(self, force=False):
        result_dict = self.clique_comparison_metric(force=force)
        total_weight = sum([result_dict[x]['numOfClusters'] for x in result_dict.keys()])
        weighted_precision = 1.0*sum([result_dict[x]['precision'] * result_dict[x]['numOfClusters'] for x in result_dict.keys()]) / total_weight

        return weighted_precision

    def clique_comparision_metric_weighted_recall(self, force=False):
        result_dict = self.clique_comparison_metric(force=force)

        weighted_recall = 1.0*sum([result_dict[x]['recall']*result_dict[x]['numOfClusters'] for x in result_dict.keys()]) / sum([result_dict[x]['numOfClusters'] for x in result_dict.keys()])
        return weighted_recall


    #kdrew: calculate precision, recall and f1score for all clique sizes up to largest cluster
    def clique_comparison_metric(self, force=False):

        #kdrew: if already computed, just return the saved result
        if self.clique_comparison_metric_results != None and not force:
            return self.clique_comparison_metric_results

        #kdrew: only evaluate on proteins that are in gold standard set
        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters()]
        #print "clusters in gold standard"
        #print clusters
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

        #print "max_len: %s" % max_len
        if max_len < 2:
            raise Gold_Standard_Overlap_Exception("ERROR: Gold Standard Overlap Warning: no pairs in clusters overlap with gold standard")

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

            if self.normalize_by_combinations:
                #kdrew: for each cluster weight true positives by the size of the cluster divided by the number of cliques
                tp_curr = tp_curr * (1.0*len(clust)/misc.comb(len(clust),clique_size))

            true_positives += tp_curr

            fp_curr = len(is_positive_list) - tp_curr 
            if self.normalize_by_combinations:
                #kdrew: for each cluster weight false positives by the size of the cluster divided by the number of cliques
                fp_curr = fp_curr * (1.0*len(clust)/misc.comb(len(clust),clique_size))

            false_positives += fp_curr

        for gs_clust in self.get_gold_standard():
            is_positive_list = [ np.max(map(set(gs_group).issubset,clusters)) for gs_group in it.combinations(gs_clust, clique_size) ]
            gs_tp_curr = sum(is_positive_list)
            if self.normalize_by_combinations:
                #kdrew: for each cluster weight gs true positives by the size of the gold_standard cluster divided by the number of cliques
                gs_tp_curr = gs_tp_curr * (1.0*len(gs_clust)/misc.comb(len(gs_clust),clique_size))
            gs_true_positives += gs_tp_curr

            fn_curr = sum(np.logical_not(is_positive_list))
            if self.normalize_by_combinations:
                #kdrew: for each cluster weight false negatives by the size of the gold_standard cluster divided by the number of cliques
                fn_curr = fn_curr * (1.0*len(gs_clust)/misc.comb(len(gs_clust),clique_size))
            false_negatives += fn_curr

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
        random_clique_weights = dict()
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
            #kdrew: for each random cluster clique, weight by the size of the random cluster divided by the number of cliques
            random_clique_weights[frozenset(shuffled_l[:clique_size])] = (1.0*len(clust)/misc.comb(len(clust),clique_size)) 

        for random_clique in random_cliques:
            #kdrew: the max just finds if any are true, could probably use .any here
            if np.max(map(random_clique.issubset,self.get_gold_standard())):
                if self.normalize_by_combinations:
                    true_positives += 1 * random_clique_weights[random_clique]
                else:
                    true_positives += 1 
            else:
                if len(self.get_exclusion_complexes()) > 0 and np.max(map(random_clique.issubset,self.get_exclusion_complexes())):
                    #print "excluded clique"
                    continue
                else:
                    #print "false_positive: %s" % (' '.join(random_clique),)
                    if self.normalize_by_combinations:
                        false_positives += 1 * random_clique_weights[random_clique]
                    else:
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
        random_gs_clique_weights = dict()
        for s in xrange(self.samples):

            #kdrew: get a random cluster
            gs_clust = gs_clusters[gs_wrg()]

            shuffled_l = rand.permutation(list(gs_clust))
            #print "sampled: %s" % (shuffled_l[:clique_size],)

            random_gs_cliques.add(frozenset(shuffled_l[:clique_size]))
            random_gs_clique_weights[frozenset(shuffled_l[:clique_size])] = (1.0*len(gs_clust)/misc.comb(len(gs_clust),clique_size)) 

        for random_clique in random_gs_cliques:

            #if np.max(map( set( shuffled_l[:clique_size] ).issubset, clusters )):
            if np.max(map( random_clique.issubset, clusters )):
                if self.normalize_by_combinations:
                    gs_true_positives += 1 * random_gs_clique_weights[random_clique]
                else:
                    gs_true_positives += 1 
            else:
                if self.normalize_by_combinations:
                    false_negatives += 1 * random_gs_clique_weights[random_clique]
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
    parser.add_argument("--remove_non_gold_standard_proteins", action="store_true", dest="remove_non_gold_standard_proteins", required=False, default=True,
                                            help="Flag to remove proteins from clusters that are not in the gold standard, default=True")
    parser.add_argument("--normalize_by_combinations", action="store_true", dest="normalize_by_combinations", required=False, default=True,
                                            help="Normalize clique precision recall by the number of combinations for each cluster, default=True")
    parser.add_argument("--excluded_complexes", action="store", dest="excluded_complexes", required=False, default=None,
                                            help="Filename of benchmark complexes to be excluded from false positive calculation, default=None")
    parser.add_argument("--pseudocount", action="store", type=float, dest="pseudocount", required=False, default=0.00001,
                                            help="Set pseudocount for clique sampling, default=0.00001")

    args = parser.parse_args()

    gold_standard_complexes = []
    gold_file = open(args.gold_standard_filename,"rb")
    for line in gold_file.readlines():
        gold_standard_complexes.append(line.split())
    gold_file.close()

    predicted_clusters = []
    clpred_f = open(args.cluster_filename,"rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())
    clpred_f.close()

    excluded_complexes = []
    if args.excluded_complexes != None:
        exclude_file = open(args.excluded_complexes,"rb")
        for line in exclude_file.readlines():
            excluded_complexes.append(line.split())

        exclude_file.close()

    cplx_compare = ComplexComparison(gold_standard_complexes, predicted_clusters, remove_non_gold_standard_proteins=args.remove_non_gold_standard_proteins, exclusion_complexes=excluded_complexes, normalize_by_combinations=args.normalize_by_combinations, pseudocount=args.pseudocount)
    #print "Sensitivity: %s" % cplx_compare.sensitivity()
    #print "PPV: %s" % cplx_compare.ppv()
    #print "ACC: %s" % cplx_compare.acc()
    #print "MMR: %s" % cplx_compare.mmr()
    #print "PWMMR: %s" % cplx_compare.pwmmr()
    #print "MMR_PWMMR_hmean: %s" % cplx_compare.mmr_pwmmr_hmean()
    #print "Precision measure: %s" % cplx_compare.precision_measure()
    #print "Recall measure: %s" % cplx_compare.recall_measure()
    #print "Precision Recall product: %s" % cplx_compare.precision_recall_product()
    #ccmm = cplx_compare.clique_comparison_metric_mean()
    #print "Clique Precision Mean: %s Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
    #ccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
    #print "Clique Weighted Precision Mean: %s Weighted Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])
    #print "Clique Weighted hmean (F-weighted K-Clique): %s" % (hmean([ccmm['precision_mean'],ccmm['recall_mean']]))

    sensitivity = cplx_compare.sensitivity()
    ppv = cplx_compare.ppv()
    acc = cplx_compare.acc()
    mmr = cplx_compare.mmr()
    pwmmr = cplx_compare.pwmmr()
    pwmmr_hmean = cplx_compare.mmr_pwmmr_hmean()
    precision_measure = cplx_compare.precision_measure()
    recall_measure = cplx_compare.recall_measure()
    precision_recall_product = cplx_compare.precision_recall_product()
    ccmm = cplx_compare.clique_comparison_metric_mean()
    clique_pr_mean = ccmm['precision_mean']
    clique_re_mean = ccmm['recall_mean']
    clique_f1grand = cplx_compare.clique_comparison_metric_grandf1score(mean_func=np.mean)
    clique_weighted_precision = cplx_compare.clique_comparision_metric_weighted_precision()
    clique_weighted_recall = cplx_compare.clique_comparision_metric_weighted_recall()
    wccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
    clique_weighted_pr_mean = wccmm['precision_mean']
    clique_weighted_re_mean = wccmm['recall_mean']
    clique_weighted_hmean = hmean([wccmm['precision_mean'],wccmm['recall_mean']])
    print "Sensitivity\tPPV\tACC\tMMR\tPWMMR\tMMR_PWMMR_hmean\tPrecision measure\tRecall measure\tPrecision Recall product\tClique Precision Mean\tRecall Mean\tF-Grand K-Clique\tClique Weighted Precision Mean\tWeighted Recall Mean\tClique Weighted hmean (F-weighted K-Clique), Clique Weighted Precision, Clique Weighted Recall\n"
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sensitivity, ppv, acc, mmr, pwmmr, pwmmr_hmean, precision_measure, recall_measure, precision_recall_product, clique_pr_mean, clique_re_mean, clique_f1grand, clique_weighted_pr_mean, clique_weighted_re_mean, clique_weighted_hmean, clique_weighted_precision, clique_weighted_recall) 


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

class Gold_Standard_Overlap_Exception(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

if __name__ == "__main__":
        main()

