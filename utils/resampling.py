import numpy as np

class FeatureResampling:

    def __init__(self): pass

    def _poisson_noise(self,df):
        '''Return dataframe with poisson noise added'''
        ## bjl: Not sure why we just add, should we randomly subtract as well?? ##
        return df + np.random.poisson(df)
        
    def _bootstrap(self,df):
        '''Sample columns with replacement, returning df of same shape'''
        return df.sample(df.shape[1], replace=True, axis=1)
