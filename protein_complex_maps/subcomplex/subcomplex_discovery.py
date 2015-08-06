
import logging
import numpy as np
import itertools as it
import multiprocessing as mp

import argparse
import pickle
import pandas as pd

import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.protein_util as pu

pd.options.display.max_columns = 50

logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

    parser = argparse.ArgumentParser(description="Discover subcomplexes in fractionation data")
    parser.add_argument("--input_msds_pickles", action="store", nargs='+', dest="msds_filenames", required=True,
                                help="Filenames of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, default=[],
                                help="Protein ids in which to anaylze")
    parser.add_argument("--protein_file", action="store", dest="protein_file", required=False, default=None,
                                help="File of protein ids in which to anaylze, each line has separate cluster of proteins")
    parser.add_argument("--subcomplex_result_pickle", action="store", dest="subcomplex_result_pickle", required=True,
                                help="Output file for which to pickle subcomplex results")
    parser.add_argument("--threshold", action="store", dest="threshold", type=float, required=False, default=0.0,
                                help="Threshold for which consider protein present or absent in fraction, default=0.0")
    parser.add_argument("--merge_threshold", action="store", dest="merge_threshold", type=float, required=False, default=0.5,
                                help="Threshold for which consider merging two clusters by jaccard index, default=0.5")
    parser.add_argument("--filter_threshold", action="store", dest="filter_threshold", type=float, required=False, default=0.5,
                                help="Threshold for which consider filter two clusters by jaccard index (only include highest zscore cluster), default=0.5")
    parser.add_argument("--output_file", action="store", dest="output_file", required=False, default=None,
                                help="Filename for resulting subcomplex proteins, sorted by zscore_all DESC")

    args = parser.parse_args()


    #scorefxn = su.multiple_dot
    #scorefxn = su.sum_matrix

    cluster_list = []
    if len(args.proteins) > 0:
        cluster_list.append(args.proteins)

    elif args.protein_file != None:
        f = open(args.protein_file,"rb")
        for line in f.readlines():
            cluster_list.append(line.split())
        f.close()

    else:
        print "Error: need to set args.protein_file or args.proteins"
        return -1


    filtered_subcomplexes = dict()
    #kdrew: go through all of the clusters in the input and calculate subcomplexes
    #for cluster_proteins in cluster_list:
    multiproc_input = [(cluster_proteins,args) for cluster_proteins in cluster_list]
    p = mp.Pool(12)
    subcomplex_results = p.map(multiproc_complex_discovery_helper, multiproc_input)

    for fs in subcomplex_results:
        filtered_subcomplexes.update(fs)

        for result in sorted(fs.values()):
            #uniprot_ids = [cd_obj.rev_protein_map[prot_id] for prot_id in result.proteins]
            #genenames = [genename_map[prot_id] for prot_id in uniprot_ids]
            #print genenames
            print result
            print "\n"

    if args.output_file != None:
        ofile = open(args.output_file,"wb")

        for result in sorted(filtered_subcomplexes.values(), reverse=True):
            ofile.write(' '.join(result.proteins))
            ofile.write("\n")

        ofile.close()

    pickle.dump( sorted(filtered_subcomplexes.values(), reverse=True), open(args.subcomplex_result_pickle,"wb"))


def multiproc_complex_discovery_helper(parameter_tuple):
    cluster_proteins = parameter_tuple[0]
    args = parameter_tuple[1]

    print cluster_proteins
    genename_map = pu.get_genenames_uniprot( cluster_proteins )

    cd_obj = ComplexDiscovery( args.msds_filenames, cluster_proteins, threshold=args.threshold, scorefxn=su.sum_matrix)
    cd_obj.discover_complexes()
    cd_obj.merge_subcomplexes(args.merge_threshold)
    cd_obj.calc_zscore()
    fs = cd_obj.filter_subcomplexes(args.filter_threshold, remove_singletons=True)

    return fs

class SubcomplexResult(object):
    def __init__(self, proteins=[], fractions=[]):
        self.proteins = frozenset(proteins)
        self.fractions = set(fractions)
        self.zscore_all = None
        self.zscore_cols = None
        self.zscore_rows = None

    def __repr__(self):
        return "proteins: %s\nfractions: %s\nzscore_all: %s\nzscore_cols: %s\nzscore_rows: %s" % (self.proteins, self.fractions, self.zscore_all, self.zscore_cols, self.zscore_rows)

    def __cmp__(self, obj):
        if obj.zscore_all == None or self.zscore_all > obj.zscore_all:
            return 1
        elif self.zscore_all == None or self.zscore_all < obj.zscore_all:
            return -1
        else:
            return 0


    #def add_proteins(prots):
    #    self.proteins = self.proteins.union(prots)
    def add_fractions(self, fracs):
        self.fractions = self.fractions.union(fracs)

class ComplexDiscovery(object):

    def __init__(self, msds_filenames, proteins, threshold=0.0, scorefxn=su.sum_matrix, normalize_flag=True):
        self.msds_filenames = msds_filenames
        self.input_proteins = proteins
        self.threshold = threshold
        self.proteins = []
        self.protein_map = dict()
        self.rev_protein_map = dict()
        self.normalize_flag = normalize_flag
        self.scorefxn = scorefxn

        self.results_list = []
        self.df_dict = dict()
        self.df_index_dict = dict()

        self.subcomplex_results = dict()

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
                    #kdrew: do the reverse mapping
                    self.rev_protein_map[name_list[id_dict[prot_id]]] = prot_id
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


            if self.normalize_flag:
                #kdrew: normalize by proteins (sum of all protein vectors = 1.0, i.e. probability)
                df = df.div( df.sum(axis=0), axis=1 )
                #print df.sum(axis=0)

            df = df.fillna(0.0)

            #kdrew: store dataframe
            self.df_dict[ msds_filename ] = df

    #kdrew: function to find combinations of proteins for each given fraction
    def discover_complexes(self,):

        for i, df_key in enumerate(self.df_dict):
            df = self.df_dict[df_key]
            #print df

            #df_complex = df[list(self.proteins)]
            #print df_complex

            for fraction in df.index:
                #print fraction
                df_frac = df.ix[fraction,list(self.proteins)]
                df_frac = df_frac[df_frac > self.threshold]
                sorted_prot_ids = sorted(df_frac.index.tolist())
                #print sorted_prot_ids
                #print df_frac
                mean_val = df_frac.mean()
                #print "mean_val: %s" % (mean_val,)
                try:
                    self.subcomplex_results[frozenset(sorted_prot_ids)].add_fractions([fraction]) 
                except KeyError:
                    self.subcomplex_results[frozenset(sorted_prot_ids)] = SubcomplexResult(proteins=sorted_prot_ids, fractions=[fraction])




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
    
    def filter_subcomplexes(self, filter_threshold=0.5, remove_singletons=False):
        filtered_subcomplexes = dict()

        #kdrew: sort by zscore_all, highest first
        for prospect_subcomplex in sorted(self.subcomplex_results.values(), reverse=True):
            #kdrew: lousy name for a flag, is true if not similar to anything else
            if remove_singletons and len(prospect_subcomplex.proteins) < 2:
                continue
            else:
                not_similar2filtered = True
            for prots in filtered_subcomplexes:
                if jaccard_index( prots, prospect_subcomplex.proteins ) > filter_threshold:
                    #kdrew: found a subcomplex that is similar to current one, do not include in final results
                    not_similar2filtered = False
                    break

            if not_similar2filtered:
                filtered_subcomplexes[frozenset(prospect_subcomplex.proteins)] = prospect_subcomplex


        return filtered_subcomplexes


    #kdrew: there are a few ways to merge complexes, merge their intersection, merge their union, merge based on probability distribution of jaccard index 
    #kdrew: this should also be implemented in the bicluster framework and allowed to sample in the monte carlo fashion
    def merge_subcomplexes(self, merge_threshold=0.5):

        for subcomplex_i in self.subcomplex_results.values():
            proteins_i = subcomplex_i.proteins
            if len(proteins_i) == 0:
                continue
            for subcomplex_j in self.subcomplex_results.values():
                proteins_j = subcomplex_j.proteins
                #print "proteins_i: %s" % (proteins_i,)
                #print "proteins_j: %s" % (proteins_j,)
                jindex = jaccard_index(proteins_i, proteins_j)
                #print "jindex: %s" % jindex
                if merge_threshold < jindex:
                    ##kdrew: merging by union
                    #sorted_prot_ids = sorted(list(proteins_i.union(proteins_j)))
                    #try:
                    #    subcomplex_results[frozenset(sorted_prot_ids)] = list(set(subcomplex_results[frozenset(sorted_prot_ids)] + i[1] + j[1]))
                    #except KeyError:
                    #    subcomplex_results[frozenset(sorted_prot_ids)] = list(set(i[1] + j[1]))

                    #kdrew: merging by intersection
                    sorted_prot_ids = sorted(list(set(proteins_i).intersection(set(proteins_j))))
                    i_fracs = self.subcomplex_results[proteins_i].fractions
                    j_fracs = self.subcomplex_results[proteins_j].fractions
                    try:
                        #subcomplex_results[frozenset(sorted_prot_ids)] = list(set(subcomplex_results[frozenset(sorted_prot_ids)] + i[1] + j[1]))

                        self.subcomplex_results[frozenset(sorted_prot_ids)].add_fractions(i_fracs.union(j_fracs))

                    except KeyError:
                        #subcomplex_results[frozenset(sorted_prot_ids)] = list(set(i[1] + j[1]))

                        self.subcomplex_results[frozenset(sorted_prot_ids)] = SubcomplexResult(proteins=sorted_prot_ids, fractions=i_fracs.union(j_fracs) ) 

        #return subcomplex_results


    def calc_zscore(self,):

        for result in self.subcomplex_results.values():
            proteins = result.proteins
            fractions = result.fractions
            #uniprot_ids = [self.rev_protein_map[prot_id] for prot_id in proteins]
            #genenames = [genename_map[prot_id] for prot_id in uniprot_ids]
            #print "proteins: %s" % (genenames, )
            #print "fractions: %s" % (fractions, )
            #print "# of fracs: %s" % (len(fractions), )
            #print "\n"
            #
            #print "self.subcomplex: %s" % (self.subcomplex_results[frozenset(proteins)],)
            #print "\n"

            #kdrew: only testing a single df so this probably needs to be corrected if passing in multiple msds on the commandline
            for df_key in self.df_dict:
                df = self.df_dict[df_key]
                df = df[self.proteins]

                rsscore_obj = rsu.RandomSamplingScore(df.as_matrix(), self.scorefxn, sample_module=np.random)
                bicluster1 = bc.BiclusterDF(rows=fractions, cols=proteins, random_module=np.random, data_frame=df)
                #print "bicluster.rows(): %s" % bicluster1.rows()
                try:
                    zscore_all = rsscore_obj.zscore_all(bicluster1.get_submatrix(df))
                    #print "zscore: %s" % (zscore_all)

                    #kdrew: proteins are fixed, shuffle fractions
                    zscore_cols = rsscore_obj.zscore_columns(df.as_matrix(), bicluster1)
                    #print "cols zscore: %s" % (zscore_cols)

                    #kdrew: fractions are fixed, shuffle proteins
                    zscore_rows = rsscore_obj.zscore_rows(df.as_matrix(), bicluster1)
                    #print "rows zscore: %s" % (zscore_rows)

                    self.subcomplex_results[frozenset(proteins)].zscore_all = zscore_all
                    self.subcomplex_results[frozenset(proteins)].zscore_cols = zscore_cols
                    self.subcomplex_results[frozenset(proteins)].zscore_rows = zscore_rows

                except rsu.EmptyColumnError, rsu.EmptyRowError:
                    continue

def jaccard_index(x, y):
    sx = set(x)
    sy = set(y)
    return 1.0 * len(sx.intersection(sy)) / len(sx.union(sy))

#kdrew: annoying kludge to get object oriented (or any multiparameter) functions to work with multiprocessing
#kdrew: outside of class, unpacks arguments including 'self' and calls class function with repackaged args
def multiprocess_helper(arg, **kwargs):
    self = arg[0][0]
    df_key = arg[0][1]
    protein_list = arg[1]
    return self.find_fractions(df_key, protein_list)

if __name__ == "__main__":
    main()

