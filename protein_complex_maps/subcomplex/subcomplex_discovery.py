
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
    parser.add_argument("--subcomplex_result_pickle", action="store", dest="subcomplex_result_pickle", required=True,
                                help="Output file for which to pickle subcomplex results")
    parser.add_argument("--maxr", action="store", type=int, dest="maxr", required=False, default=10,
                                help="When calculating protein sets using combinations what is the max choose r to evaluate, will find all sets from pairs to size maxr, default=10")
    parser.add_argument("--threshold", action="store", dest="threshold", type=float, required=False, default=0.0,
                                help="Threshold for which consider protein present or absent in fraction, default=0.0")

    args = parser.parse_args()

    cd_obj = ComplexDiscovery( args.msds_filenames, args.proteins, threshold=args.threshold)
    subcomplex_results = cd_obj.discover_complexes(maxr=args.maxr)

    pickle.dump(subcomplex_results, open(args.subcomplex_result_pickle,"wb"))

    for i in sorted(subcomplex_results.items(), key = lambda x: len(x[1][0][1])):
        print i, len(i[1][0][1])


class ComplexDiscovery(object):

    def __init__(self, msds_filenames, proteins, threshold=0.0, normalize_flag=True):
        self.msds_filenames = msds_filenames
        self.input_proteins = proteins
        self.threshold = threshold
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
                #kdrew: normalize by proteins (sum of all protein vectors = 1.0, i.e. probability)
                df = df.div( df.sum(axis=0), axis=1 )
                #print df.sum(axis=0)

            #kdrew: threshold out fraction vectors that do not have enough of the input proteins
            #print self.proteins
            #print df
            #bool_vector = df[self.proteins] > 0.0
            #sum_vector = bool_vector.sum(axis=1)
            #percent_vector = 1.0*(sum_vector)/len(self.proteins)
            #percent_vector = 1.0*(df[self.proteins] > 0.0).sum(axis=1)/len(self.proteins)
            #df = df[percent_vector > self.protein_percent_threshold]


            #kdrew: store list of fractions that pass threshold
            #self.df_index_dict[msds_filename] = list(df[percent_vector > self.protein_percent_threshold].index)

            #kdrew: store dataframe
            self.df_dict[ msds_filename ] = df

    #kdrew: multiprocessor function to find fractions for all combinations of proteins sets
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

                results = pool.map(multiprocess_helper, it.izip_longest([],it.combinations(self.proteins,r),fillvalue=(self,df_key)))

                for res in results:
                    if len(res[2]) > 0:
                        protein_list = res[0]
                        df_key = res[1]
                        pl_only_fractions = res[2]
                        print "protein_list: %s msds: %s fractions: %s" % (protein_list, df_key, pl_only_fractions)
                        print df.ix[pl_only_fractions,list(protein_list)]
                        try:
                            subcomplex_dict[protein_list].append((df_key, pl_only_fractions))
                        except KeyError:
                            subcomplex_dict[protein_list] = [(df_key,pl_only_fractions)]

        return subcomplex_dict


    #kdrew: function to find fractions that only have proteins in protein_list above some threshold (default 0.0) and no other protein above threshold
    def find_fractions(self, df_key, protein_list):
        df = self.df_dict[df_key]
        #kdrew: protein_list is a list of proteins that are potentially members of a subcomplex
        df_pl = df[list(protein_list)]
        df_pl_greater0 = df_pl[df_pl > self.threshold].dropna()
        pl_fractions = df_pl_greater0.index

        out_list = list(set(self.proteins) - set(protein_list))

        df_out = df[out_list]
        df_out_equal0 = df_out[df_out <= self.threshold].dropna()
        out_fractions = df_out_equal0.index

        pl_only_fractions = list(set(pl_fractions).intersection(set(out_fractions)))

        return (protein_list, df_key, pl_only_fractions)

#kdrew: annoying kludge to get object oriented (or any multiparameter) functions to work with multiprocessing
#kdrew: outside of class, unpacks arguments including 'self' and calls class function with repackaged args
def multiprocess_helper(arg, **kwargs):
    self = arg[0][0]
    df_key = arg[0][1]
    protein_list = arg[1]
    return self.find_fractions(df_key, protein_list)

if __name__ == "__main__":
    main()

