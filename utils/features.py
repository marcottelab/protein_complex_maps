import numpy as np
import scipy.stats as stats
import scipy.spatial.distance as dist
import math

def js_pairs(P,Q,distance=False):
    '''Compute Jensen-Shannon distance metric for two frequency vectors. No zeros allowed!'''
    kl_distance = lambda x,y: np.sum(P * np.log2(P/Q)) # this was fastest, beating a map/zip approach and scipy.entropy
    A = np.mean([P,Q],axis=0) # Calculate a mean of the two vectors
    js_diverg = .5 * (kl_distance(P,A)) + .5 * (kl_distance(Q,A))
    if distance: # return the distance (a metric), the square root of the divergence
        return math.sqrt(js_diverg)
    return js_diverg # Better to return the divergence, which is bounded (0,1)

class FeatureFunctions:

    '''All return a NxN matrix of features'''

    def __init__(self): pass
        
    def _kullback_leibler(self,df):
        '''Return matrix of Kullback_Leibler divergences.
        It's better to use Jensen-Shannon because it's symmetric'''
        return dist.squareform( dist.pdist(df, lambda x,y: np.sum(P * np.log2(P/Q))) )
    
    def _jensen_shannon(self,df):
        '''Return matrix of Jensen-Shannon divergences'''
        return dist.squareform( dist.pdist(df, lambda x,y: js_pairs(x,y)) )
        
    def _pearsonR(self,df):
        '''Return pearson correlation matrix'''
        return np.nan_to_num( np.corrcoef(df),0 )
        
    def _spearmanR(self,df):
        '''Return spearman ranked correlation coeffienct matrix'''
        rho,pval = stats.spearmanr(df.T)
        rho = np.nan_to_num(rho,0)
        return rho
        
    def _spearmanR_weighted(self,df):
        '''Return spearman ranked correlation coeffienct matrix, weighted by the p-value'''
        rho,pval = stats.spearmanr(df.T)
        rho = np.nan_to_num(rho,0)
        pval = np.nan_to_num(pval,0)
        return rho * (1 - pval)
        
    def _euclidean(self,df):
        '''Return euclidean distance'''
        return dist.squareform(dist.pdist(df,'euclidean'))
        
    def _covariance(self,df):
        '''Return the covariance between two input arrays'''
        return np.nan_to_num( np.cov(df),0 )