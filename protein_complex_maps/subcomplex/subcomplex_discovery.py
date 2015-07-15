
import logging
import numpy as np
import itertools as it
import multiprocessing as mp

import argparse
import pickle
import pandas as pd

pd.options.display.max_columns = 50

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

    parser = argparse.ArgumentParser(description="Discover subcomplexes in fractionation data")
    parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True,
                                help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True,
                                help="Protein ids in which to anaylze")
    #parser.add_argument("--protein_percent_threshold", action="store", dest="protein_percent_threshold", type=float, required=False, default=0.9,
    #                            help="Percentage of proteins that need to be present in input fractions (expressed as decimal), default 0.9 (ie 90%)")

    args = parser.parse_args()

    cd_obj = ComplexDiscovery( args.msds_filenames, args.proteins)
    cd_obj.discover_complexes(maxr=4)


class ComplexDiscovery(object):

    def __init__(self, msds_filenames, proteins, normalize_flag=True):
        self.msds_filenames = msds_filenames
        self.input_proteins = proteins
        self.proteins = []
        self.protein_map = dict()
        self.normalize_flag = normalize_flag

        self.results_list = []
        self.df_dict = dict()
        self.df_index_dict = dict()

        #kdrew: initialize
        for prot_id in self.input_proteins:
            self.protein_map[prot_id] = prot_id

        self.create_data_frames()


    def create_data_frames(self,):

        msds_dict= dict()
        for msds_filename in self.msds_filenames:
            msds = pickle.load( open( msds_filename, "rb" ) )
            msds_dict[msds_filename] = msds

            #kdrew: map input protein ids to protein ids used for indexing
            #kdrew: id_dict is the mapping from protein_id -> index
            id_dict = msds.get_id_dict()
            #kdrew: name list will map id -> protein_id used for indexing
            name_list = msds.get_name_list()
            #kdrew: this will fill protein_map with the id used in msds
            for prot_id in self.input_proteins:
                try:
                    self.protein_map[prot_id] = name_list[id_dict[prot_id]]
                except KeyError:
                    continue

        self.proteins = self.protein_map.values()
        print "self.proteins"
        print self.proteins

        #kdrew: store data in data frames
        for msds_filename in msds_dict.keys():
            msds = msds_dict[msds_filename]
            df = msds.get_data_frame()

            #print "filename: %s" % msds_filename

            missing_proteins = [ self.protein_map[i] for i in self.input_proteins if self.protein_map[i] not in df.columns ]

            print "final mapped proteins: %s" % self.proteins
            print "final missing proteins: %s" % missing_proteins

            #kdrew: add in 0.0 rows for missing entries that are still unmapped
            for missing_prot in missing_proteins:
                df[missing_prot] = np.repeat( 0.0, len(df.index) )


            ##kdrew: normalize by fractions (max of all fraction vectors = 1.0)
            #df = df.div( df.max(axis=1), axis=0 )

            if self.normalize_flag:
                #kdrew: normalize by fractions (sum of all fraction vectors = 1.0, i.e. probability)
                df = df.div( df.sum(axis=1), axis=0 )

            #kdrew: threshold out fraction vectors that do not have enough of the input proteins
            #print self.proteins
            #print df
            bool_vector = df[self.proteins] > 0.0
            sum_vector = bool_vector.sum(axis=1)
            percent_vector = 1.0*(sum_vector)/len(self.proteins)
            #percent_vector = 1.0*(df[self.proteins] > 0.0).sum(axis=1)/len(self.proteins)
            #df = df[percent_vector > self.protein_percent_threshold]


            #kdrew: store list of fractions that pass threshold
            #self.df_index_dict[msds_filename] = list(df[percent_vector > self.protein_percent_threshold].index)

            #kdrew: store dataframe
            self.df_dict[ msds_filename ] = df

    def discover_complexes(self, maxr=None):

        subcomplex_dict = dict()

        pool = mp.Pool()

        for df_key in self.df_dict:
            df = self.df_dict[df_key]
            #print df

            if maxr == None or maxr > len(self.proteins):
                maxr = len(self.proteins)

            #kdrew: for different possible numbers of sequential experiments, focusing on 1 exp or combinations of 2 or more ...
            for r in range(2, maxr+1):
                #print "r: %s" % r

                for protein_list in it.combinations( self.proteins, r ):
                    #kdrew: protein_list is a list of proteins that are potentially members of a subcomplex
                    #print protein_list
                    df_pl = df[list(protein_list)]
                    df_pl_greater0 = df_pl[df_pl > 0.0].dropna()
                    #print df_pl_greater0.head()
                    pl_fractions = df_pl_greater0.index
                    #print pl_fractions

                    out_list = list(set(self.proteins) - set(protein_list))

                    #print out_list
                    df_out = df[out_list]
                    df_out_equal0 = df_out[df_out == 0.0].dropna()
                    #print df_out_equal0.head()
                    out_fractions = df_out_equal0.index
                    #print out_fractions

                    pl_only_fractions = list(set(pl_fractions).intersection(set(out_fractions)))
                    #print "pl_only_fractions"
                    #print  pl_only_fractions

                    if len(pl_only_fractions) > 0:
                        print "protein_list: %s msds: %s fractions: %s" % (protein_list, df_key, pl_only_fractions)
                        try:
                            subcomplex_dict[protein_list].append((df_key, pl_only_fractions))
                        except KeyError:
                            subcomplex_dict[protein_list] = [(df_key,pl_only_fractions)]


if __name__ == "__main__":
    main()

