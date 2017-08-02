import numpy as np
import scipy.stats as stats
import scipy.spatial.distances as dist
import math

class FeatureFunctions:

    def __init__(self): pass

    def _kullback_leibler(self,P,Q):
        '''Compute KL divergence for two frequency vectors. No zeros allowed!'''
        
        # This needs to be done somewhere, but I think upstream is better
        #assert 0 not in P, "No zeros allowed in input vector"
        #assert 0 not in Q, "No zeros allowed in input vector"
        
        # Convert to probability dists
        P = np.array(P,float) + 1 # add 1 to all to prevent zeros
        Q = np.array(Q,float) + 1
        P = P/P.sum()
        Q = Q/Q.sum()
        
        return np.sum(P * np.log2(P/Q)) # this was fastest, beating a map/zip approach and scipy.entropy
        
    def _jensen_shannon(self,P,Q,distance=False):
        '''Compute Jensen-Shannon distance metric for two frequency vectors. No zeros allowed!'''
        A = np.mean([P,Q],axis=0) # Calculate a mean of the two vectors
        js_diverg = .5 * (_kullback_leibler(P,A)) + .5 * (_kullback_leibler(Q,A))
        if distance: # return the distance (a metric), the square root of the divergence
            return math.sqrt(js_diverg)
        return js_diverg # Better to return the divergence, which is bounded (0,1)
        
    def _pearsonR(self,P,Q):
        '''Return pearson correlation coeffienct'''
        rho,pval = stats.pearsonr(P,Q)
        return rho
        
    def _pearsonR_weighted(self,P,Q):
        '''Return pearson correlation coeffienct, weighted by the p-value'''
        rho,pval = stats.pearsonr(P,Q)
        return rho * (1 - pval)
        
    def _spearmanR(self,P,Q,weighted=True):
        '''Return spearman ranked correlation coeffienct'''
        rho,pval = stats.spearmanr(P,Q)
        return rho
        
    def _spearmanR_weighted(self,P,Q,weighted=True):
        '''Return spearman ranked correlation coeffienct, weighted by the p-value'''
        rho,pval = stats.spearmanr(P,Q)
        return rho * (1 - pval)
        
    def _euclidean(self,P,Q):
        '''Return euclidean distance'''
        return dist.euclidean(P,Q)
        
    def _covariance(self,P,Q):
        '''Return the covariance between two input arrays'''
        return np.cov(P,Q)[0,1]