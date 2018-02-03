#! /usr/bin/env python

import numpy as np
import pandas as pd
import itertools as it

from functions import features, resampling


class Elut():

    '''
    Class to read and inspect a fractionaton masss-spec experiment
    
    Infile must be wide format with the first column holding protein or orthogroup ids.
    '''
    
    def __init__(self,data=None):
        
        ### Could just extend the DataFrame object by inheriting it, doesn't feel right, though ###
        # self._df = pd.DataFrame.__init__(self,data)
        
        self.df = None
        self.info = {}
        if not data is None:
            self.df = data
            self._check_dtypes()
            self._get_info()
            
        self.is_normalized = False
        self.is_thresholded = False
            
        ### Write attributes to hold data that can be shown with self.description. Also generate
        ### values for pseudocounts that go into poisson noise and kullback-leibler and jensen-shannon
        ### features
        
    def _check_dtypes(self):
        '''Make sure loaded DataFrame is only floats'''
        try:
            dtype_set = set(self.df.dtypes)
            assert len(dtype_set) == 1
            assert dtype_set.pop() is np.dtype("float")
        except AssertionError:
            raise Exception("Multiple datatypes or non-float datatypes found in input DataFrame")
    
    def load(self,infile,format='csv'):
        '''Read in data as pandas dataframe in wide format'''
        if not self.df is None:
            raise Exception("data already loaded")
        if format == 'csv':
            self.df = pd.read_csv(infile,index_col=0).astype("float")
        elif format == 'csv':
            self.df = pd.read_table(infile,index_col=0).astype("float")
        else:
            raise Exception("<format> must be either 'csv' or 'tsv'")
        self._check_dtypes()
        self._get_info()
        
    def load_many(self,file_list):
        '''Load a list of files and join together (by row) into a combined dataframe'''
        ## TBD
        pass
            
    def normalize(self,by='row'):
        '''Normalize the loaded dataframe by row or by column sums'''
        assert not self.is_normalized, "DataFrame already normalized"
        self._prenormed_df = self.df # save df in hidden var before normalizing
        if by == 'row':
            self.df = self.df.div(self.df.sum(axis=1), axis=0)
        elif by == 'column':
            self.df = self.df/self.df.sum()
        else: raise Exception("must specify either 'row' or 'column'")
        self.is_normalized = True
        
    def threshold(self,thresh=5):
        '''Drop rows with total PSMs < thresh'''
        assert not self.is_thresholded, "DataFrame already thresholded"
        self.df["sum"] = self.row_sums
        len_before = len(self.df)
        self.df = self.df[ self.df["sum"] >= thresh ]
        self.df.drop("sum",axis=1,inplace=True)
        print "Dropped {} rows".format(len_before - len(self.df))
        self.is_thresholded = True
        
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
        
    def _get_info(self):
        '''Output simple info on the stored elution profiles'''
        if self.df is None:
            print "No loaded data"
            return
        
        self.row_sums = self.df.sum(axis=1)
        self.info["n_PSMs"] = self.row_sums.sum()
        self.info["PSMs_per_row_stats"] = self.row_sums.describe()[1:].to_dict()
        self.info["n_proteins"] = len(self.df)
        self.info["n_fractions"] = len(self.df.columns)
        

        
            
class ElutFeatures(Elut,features.FeatureFunctions,resampling.FeatureResampling):

    '''Class to extract a variety of quantitative features on each pair of
    proteins/orthogroups in an elution experiment'''
    
    ### Notes:
    ###     - Would be nice to add some meta data on bounds etc. for the features
    
    available_features = ["pearsonR",
                "spearmanR",
                "spearmanR_weighted",
                "jensen_shannon",
                "kullback_leibler",
                "euclidean",
                "covariance",
                "xcorr"]
                
    resampling_strategies = ["poisson_noise",
                            "bootstrap"]
    
    def __init__(self,data=None):
        Elut.__init__(self,data)
        
        self.analysis_count = 0
        self.analyses = []

            
    def _to_df(self,df,feature_matrix,feature_string):
        '''Turn the 1d output of pdist into a tidy DataFrame'''
        ## Return off diagonal by default - make this a flag? ##
        square_df = pd.DataFrame( feature_matrix, columns=df.index.values, index=df.index )
        square_df.sort_index(inplace=True)
        square_df.sort_index(axis=1,inplace=True)
        mask = np.triu( np.ones(square_df.shape),k=1 ).astype('bool')
        off_diag = square_df.where(mask)
        tidy_df = off_diag.stack().reset_index()
        tidy_df.columns = ["ID1", "ID2", feature_string]
        return tidy_df
    
    def _average_resamples(self,df,feature_function,resample_function,iterations):
        '''Average N resamples'''
        ## work on 1d pdist output
        return ( sum( feature_function( resample_function(df,rep=i) ) for i in range(iterations) ) / iterations )
        
    def extract_features(self,feature,resampling=None,iterations=None,threshold=None):
        '''Return a DataFrame of features'''
        
        # Assign the function to extract features
        feat = "_" + feature # because I'm hiding actual feature functions
        assert hasattr(self,feat), "{} not in available features:\n{}".format(feature,available_features)
        feat_func = getattr(self,feat) # lookup the relvant function
        
        self._analysis = {"feature": feature} # initialize or reset dictionary to hold information on which analysis was done
        
        # Assign and execute normalization functions
        ## Currently in a hacky state ##
        
        if threshold:
            self.threshold(threshold)
            try:
                assert len(self.df) > 0
            except AssertionError:
                raise Warning("No rows left after thresholding")
                return None
            self._analysis["threshold"] = "thresh"+str(threshold)
        
        if feature in ["jensen_shannon","kullback_leibler"]:
            if self.is_normalized:
                raise Warning('''Don't normalize by rows before extracting Kullback-Leibler or 
                Jensen-Shannon features. It will be done again.''')
            print "Adding pseudocounts and normalizing for JS or KL divergence"
            self.df = self.df + 1
            self.normalize(by='row')
            
        ## Don't like this part, because it could add pseudocounts twice if jensen-shannon + poisson
        if resampling == "poisson_noise": # this is how Traver's function does it, for some reason
            self.df = self.df + (1/len(self.df))
        
        # Assign resampling function and return averaged DataFrame
        if resampling:
            
            assert iterations > 1, "if resampling, must specify more than 1 iteration"
            respl = "_" + resampling
            assert hasattr(self,respl), "{} not in available resampling strategies:\n{}".format(feature,resampling_strategies)
            respl_func = getattr(self,respl)
            
            feature_matrix = self._average_resamples(self.df,feat_func,respl_func,iterations)
            
            self._analysis["resampling"] = resampling
            self._analysis["reps"] = str(iterations)+"reps"
        
        else:# If no resampling execute feature function and return DataFrame
            feature_matrix = feat_func(self.df)
        
        # Store information on what analysis was done
        a_list = []
        for a in ["feature", "resampling", "reps", "threshold"]:
            if a in self._analysis: 
                a_list.append(self._analysis[a])
        self.analyses.append( "_".join(a_list) )
        self.analysis_count += 1
        
        return self._to_df(self.df,feature_matrix,feature)
