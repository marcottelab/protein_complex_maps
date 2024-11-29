import numpy as np
import scipy.stats as stats
import scipy.spatial.distance as dist
import math

def js_pairs(P,Q,distance=False):
    '''Compute Jensen-Shannon distance metric for two frequency vectors. No zeros allowed!'''
    kl_distance = lambda x,y: np.sum(P * np.log10(P/Q)) # this was fastest, beating a map/zip approach and scipy.entropy
    A = (P + Q) / 2 # Calculate a mean of the two vectors
    js_diverg = .5 * (kl_distance(P,A)) + .5 * (kl_distance(Q,A))
    if distance: # return the distance (a metric), the square root of the divergence
        return math.sqrt(js_diverg)
    return js_diverg # Better to return the divergence, which is bounded (0,1)

def sum_difference_pairs(P,Q):
    '''Compute sum of the absolute value differences for two frequency vectors.'''
    D = P-Q
    abs_D = np.abs(D)
    sum_abs_D = np.sum(abs_D)
    return sum_abs_D

def xcorr(P,Q, zscore=True):
    '''Compute normalized cross correlation for two vectors
       If zscor==True, will calculate (X - mean(all_lags)) / stdev(all_lags)
       Otherwise, it calculates only zero-lag value by default

       np.correlate with 'full' arguments and for ex. lag = 2 returns vector with order:
             [-2lag, -1lag, 0lag, 1lag, 2lag], so the middle value is the zero lag correlation.
    '''
    if zscore ==  True:
        a = np.correlate(P, Q, "full")
        xcorrval = (a[len(a)/2] - np.mean(a)) / np.std(a)

    elif zscore == False:
        xcorrval = np.correlate(P, Q)

    return xcorrval

class FeatureFunctions:

    '''All return a NxN matrix of features'''

    def __init__(self): pass
    
    def _jensen_shannon(self,df):
        '''Return matrix of Jensen-Shannon divergences'''
        return dist.squareform( dist.pdist(df, lambda x,y: js_pairs(x,y)) )
        
    def _xcorr(self,df):
        '''Return cross correlation matrix'''
        return dist.squareform( dist.pdist(df, lambda x,y: xcorr(x,y)) )

    def _pearsonR(self,df):
        '''Return pearson correlation matrix'''
        return np.nan_to_num( np.corrcoef(df) )
        
    def _spearmanR(self,df):
        '''Return spearman ranked correlation coefficient matrix'''
        rho,pval = stats.spearmanr(df.T)
        rho = np.nan_to_num(rho)
        return rho
        
    def _spearmanR_weighted(self,df):
        '''Return spearman ranked correlation coefficient matrix, weighted by the p-value'''
        rho,pval = stats.spearmanr(df.T)
        rho = np.nan_to_num(rho)
        pval = np.nan_to_num(pval)
        return rho * (1 - pval)
        
    def _euclidean(self,df):
        '''Return euclidean distance'''
        return dist.squareform(dist.pdist(df,'euclidean'))
        
    def _covariance(self,df):
        '''Return the covariance between two input arrays'''
        return np.nan_to_num( np.cov(df) )
        
    def _canberra(self,df):
        '''Return canberra distance matrix'''
        return dist.squareform(dist.pdist(df,'canberra'))
        
    def _braycurtis(self,df):
        '''Return canberra distance matrix'''
        return dist.squareform(dist.pdist(df,'braycurtis'))
        
    def _invbraycurtis(self,df):
        '''Return canberra distance matrix'''
        return 1. - dist.squareform(dist.pdist(df,'braycurtis'))
        
    def _cosine(self,df):
        '''Return the cosine distance matrix'''
        return dist.squareform(dist.pdist(df,'cosine'))

    def _sum_difference(self,df):
        '''Return matrix of sums of absolute value differences between vectors'''
        return dist.squareform( dist.pdist(df, lambda x,y: sum_difference_pairs(x,y)) )



