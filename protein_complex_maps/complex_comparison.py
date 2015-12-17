    
import numpy as np
import pandas as pd
import pickle as p
import argparse
import itertools as it


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
        df = self.get_na_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        #print df
        max_df = df.max(axis=1)
        #print max_df
        sum_max_df = max_df.sum()
        #print sum_max_df
        num_of_gold_complexes = len(self.get_gold_standard())
        #print num_of_gold_complexes

        mmr = sum_max_df / num_of_gold_complexes
        #print mmr
        return mmr



    #kdrew: method from Brohee and Helden BMC Bioinformatics, http://www.biomedcentral.com/content/pdf/1471-2105-7-488.pdf
    #kdrew: intersection/size_of_gold = percent coverage of given complex
    #kdrew: find max
    #kdrew: size_of_gold * max_percent_coverage  / size_of_gold
    def sensitivity(self, topN=None):
        df = self.get_intersection_table()
        #kdrew: if we want to only do a range of clusters (topN), we probably want to slice sn_df here
        if topN != None:
            df = df.ix[:,:topN]
        sn_df = df.divide(map(len,self.get_gold_standard()), axis=0)
        #print sn_df.max()
        sn_co = sn_df.max(axis=1)
        #print sn_co.max()
        goldstd_sizes = map(len,self.get_gold_standard())
        Sn = 1.0*sum(sn_co * goldstd_sizes) / sum(goldstd_sizes)
        return Sn
        
    #kdrew: method from Brohee and Helden BMC Bioinformatics, http://www.biomedcentral.com/content/pdf/1471-2105-7-488.pdf
    def ppv(self, topN=None):
        df = self.get_intersection_table()
        #kdrew: if want to do a range of clusters (topN) we probably want to slice df here
        if topN != None:
            df = df.ix[:,:topN]
        ppv_df = df.divide(1.0*df.sum())
        ppv_df = ppv_df.fillna(0.0)
        PPV_cl = ppv_df.max()
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
        
        print "mean precision %s mean recall %s" % (np.mean(precision_list), np.mean(recall_list))

        return {'precision_mean':np.mean(precision_list),'recall_mean':np.mean(recall_list)}

    def clique_comparison_metric(self,):
        max_len = np.max(map(len,self.get_clusters()))
        return_dict = dict()
        for size in range(2, max_len):
            result_dict = self.clique_comparison(size)
            print result_dict
            recall = 1.0*result_dict['tp'] / (result_dict['tp'] + result_dict['fn'])
            precision = 1.0*result_dict['tp'] / (result_dict['tp'] + result_dict['fp'])
            print "precision: %s recall: %s" % (precision, recall)
            return_dict[size] = {'precision':precision,'recall':recall}

        return return_dict


    #kdrew: return chunks of size n from list l
    def chunks(l,n):
        n = max(1, n)
        #kdrew: generator that slices list into n chunks, 
        #kdrew: if the last slice is smaller than n, ignore
        return [l[i:i + n] for i in xrange(0, len(l), n) if i+n< len(l)]

    #kdrew: generator to return m number of random combinations of size n from list l
    def rand_combinations(l,n,m):
        #kdrew: determine the number of shuffling needed
        i = int(np.ceil(1.0*m/(len(l)/n)))
        for ii in xrange(i):
            shuffled_l = rand.permutation(l)
            for ch in chunks(shuffled_l,n):
                yield ch


    #kdrew: calculate confusion matrix between predicted clusters and gold standard complexes for specific clique sizes
    def clique_comparison(self, clique_size):
        true_positives = 0
        gs_true_positives = 0
        false_positives = 0
        false_negatives = 0

        for clust in self.get_clusters():
            #print "clust: %s" % (clust,)
            #kdrew: only look at cluster ids that are in the gold standard
            clust_intersection = clust & self.get_gold_standard_proteins()
            #print "clust_intersection: %s" % (clust_intersection,)
            #kdrew: if all proteins are outside of the gold standard, move on
            if len(clust_intersection) <= 0:
                continue


            is_positive_list = [ np.max(map(set(group).issubset,self.get_gold_standard())) for group in it.combinations(clust_intersection, clique_size)]
            tp_curr = sum(is_positive_list)
            true_positives += tp_curr
            #false_positives += sum(np.logical_not(is_positive_list))

            false_positives += len(is_positive_list) - tp_curr

            ##kdrew: for all combinations in cluster of size clique_size
            #for group in it.combinations(clust_intersection, clique_size):
            #    #print "group: %s" % (group,)
            #    #kdrew: check all gold standard complexes if combination is a subset  
            #    if np.max(map(set(group).issubset,self.get_gold_standard())):
            #        #print group
            #        #kdrew: if a subset of any gold standard complex, add a true positive
            #        true_positives +=1
            #    else:
            #        #kdrew: if not a subset add a false positive
            #        #kdrew: there is a guarantee here from the cluster intersection with all gold standard proteins (above), 
            #        #kdrew: that a false positive will be a set of proteins in the gold standard but not in the same complex
            #        #kdrew: this allows for extra subunits and novel clusters with no overlap with gold standard to not be counted in evaluation
            #        false_positives +=1


        for gs_clust in self.get_gold_standard():
            is_positive_list = [ np.max(map(set(gs_group).issubset,self.get_clusters())) for gs_group in it.combinations(gs_clust, clique_size) ]
            #gs_true_positives += sum(is_positive_list)
            false_negatives += sum(np.logical_not(is_positive_list))

            #for gs_group in it.combinations(gs_clust, clique_size):
            #    if np.max(map(set(gs_group).issubset,self.get_clusters())).any():
            #        #print group
            #        gs_true_positives +=1
            #    else:
            #        false_negatives +=1

        print "truepos: %s gs_truepos: %s falsepos: %s falseneg: %s" % (true_positives, gs_true_positives, false_positives, false_negatives)

        #assert true_positives == gs_true_positives

        return_dict = dict()
        return_dict['tp'] = true_positives
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

def main():

    parser = argparse.ArgumentParser(description="Compare cluster predictions to gold standard complexes")
    parser.add_argument("--cluster_predictions", action="store", dest="cluster_filename", required=True, 
                                            help="Filename of cluster predictions, format one cluster per line, ids space separated")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True, 
                                            help="Filename of gold standard complexes, format one complex per line, ids space separated")

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
    print "Precision Mean: %s Recall Mean: %s" % (ccmm['precision_mean'],ccmm['recall_mean'])


if __name__ == "__main__":
        main()

