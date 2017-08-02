#! /usr/bin/env python

import pandas as pd
import itertools as it

from utils import features, resampling


class Elut():

    '''
    Class to hold an elution profile experiment
    
    Infile must be wide format and the first column, which holds the protein/OG names, should be titled "ID"
    '''
    
    def __init__(self,data=None):
        
        ### Could just extend the DataFrame object by inheriting it, doesn't feel right, though ###
        # self._df = pd.DataFrame.__init__(self,data)
        
        self.df = None
        if not data is None:
            self.df = data
            
        ### Write attributes to hold data that can be shown with self.description. Also generate
        ### values for pseudocounts that go into poisson noise and kullback-leibler and jensen-shannon
        ### features
    
    def load(self,infile,format='csv'):
        '''Read in data as pandas dataframe in wide format'''
        if not self.df is None:
            raise Exception("data already loaded")
        if format == 'csv':
            self.df = pd.read_csv(infile) # don't make ID column index...for now
        elif format == 'csv':
            self.df = pd.read_table(infile)
        else:
            raise Exception("<format> must be either 'csv' or 'tsv'")
        assert self.df.columns[0] == "ID", "First column must be labeled 'ID'"
            
    def normalize(self,by='row'):
        '''Normalize the loaded dataframe by row or by column sums'''
        self._prenormed_df = self.df # save df in hidden var before normalizing
        if by == 'row':
            self.df = df.div(df.sum(axis=1), axis=0)
        elif by == 'column':
            self.df = df/df.sum()
        else: raise Exception("<by> must be either 'row' or 'column'")
        
    def make_tidy(self,just_return=False):
        '''Use pandas melt to reshape dataframe into tidy format, with columns "ID, "FractionID," "Total_SpecCounts".
        If <just_return> is set, it will return the tidy dataframe without changing the stored data '''
        tidy = self.df.melt(id_vars="ID",value_vars=self.df.columns[1:],var_name="FractionID",value_name="Total_SpecCounts")
        if just_return:
            return tidy
        self._wide_df = self.df # save df in hidden var before melting
        self.df = tidy
        
    def undo_tidy(self):
        '''Resets data to wide format'''
        self.df = self._wide_df
        self._wide_df = None
        
    def description(self):
        '''Output simple info on the stored elution profiles'''
        # n fractions
        # sum psms, per row, total etc.
        # n proteins
        # something about n << p problems?
        pass
        
            
class ElutFeatures(Elut,features.FeatureFunctions,resampling.FeatureResampling):

    '''Class to extract a variety of quantitative features on each pair of
    proteins/orthogroups in an elution experiment'''
    
    ### Notes:
    ###     - Would be nice to add some meta data on bounds etc. for the features
    
    def __init__(self,data=None):
        Elut.__init__(self,data)
        
        self.features_extracted = 0
        self.resampling_done = None
        
        self.available_features = ["pearsonR",
                        "pearsonR_weighted"
                        "spearmanR",
                        "spearmanR_weighted",
                        "jensen_shannon",
                        "kullback_leibler",
                        "euclidean",
                        "covariance"]
                        
        self.resampling_strategies = ["poisson_noise",
                                    "bootstrap"]
      
    
    def _feature_gen(self,df,feature_string):
        '''Return a generator of (prot1,prot2,feature)'''
        if not df.index.name == "ID": df.set_index("ID",inplace=True)
        
        # Kinda the main part
        func = getattr(self,feature_string) # Lookup function in self. All are imported from features
                                            # but hidden from users: they will be e.g. _euclidean.
                                            # So to add new features, write a function in features.py
                                            # called _my_feature or assign it to an instance of this class
                                            # but make sure it starts with an underscore (a little weird maybe?)
        
        # Now loop over all combinations of rows in the loaded data and calculate the metric
        count = 0
        for i,j in it.combinations(df.index.values,2):
            p,q = df.iloc[i].values, df.iloc[j].values # convert to arrays
            yield i,j,func(p,q)
            count += 1
            if count % 100 == 0: print count
            
    def _to_df(self,df,feature_string):
        '''Collect output of _feature_gen to a pandas DataFrame in tidy format: cols=("ID1","ID2","my_feature")'''
        rows = [row for row in self._feature_gen(df,feature_string)]
        return pd.DataFrame(rows, columns = ["ID1","ID2",feature_string)
    
    def _average_resamples(self,df,feature,resample_string,iterations):
        '''Average N resamples'''
        
        # Do i want to have these getattr calls all in extract features? maybe...
        resample_fun = getattr(self, resample_string) # Lookup resampling strategy in self
        
        
        resample_gen = ()
        return sum( self._to_df()
        
    def extract_features(self,feature,resampling=None,iterations=None):
        '''Return a DataFrame of features'''
        feat = "_" + feature # because I'm hiding actual feature functions
        assert hasattr(self,feat), "{} not in available features:\n{}".format(feature,self.available_features)
        
        if resampling:
            respl = "_" + resampling
            assert hasattr(self,respl), "{} not in available resampling strategies:\n{}".format(feature,self.resampling_strategies)
            assert iterations > 1, "if resampling, must specify more than 1 iteration"
            
            return self._average_resamples(feat,respl,iterations)
        
        else:
            self.features_extracted = [feature]
            return self._to_df(self.df, feat)
        
        
        
    #def _assign_tasks(self,features=["pearsonr"],resampling=None,iterations=None,norm=None):
    #    '''For building functions with multiple combinations of features and resampling protocols.'''
    #    
    #    ## Make sure input isn't f-ed up ##
    #    assert type(features) == list and len(features) > 0, "Must define at least one feature type: \n{}".format(self.available_features)
    #    for f in features:
    #        assert f in self.available_features, "{} not in avaiable features:\n{}".format(f,self.available_features)
    #    if resampling:
    #        assert type(resampling) == list, "Supply resampling strategies as a list"
    #        assert iterations != None, "Must specify number of iterations for resampling"
    #        for r in resampling:
    #            assert r in self.resampling_strategies, "{} not in available resampling strategies:\n{}".format(r,self.resampling_strategies)
    #            
    #        # Do something like this
    #        for feat,respl in it.product(features,resampling):
    #            self.tasks = {}
    #            self.tasks[(feat,respl): self.feature_generator(feat,respl)]
                
         