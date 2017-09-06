
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import argparse
import pickle
import numpy.random as random 


# Import nimfa library entry point for factorization
import nimfa

import logging

import protein_complex_maps.normalization_util as nu
import protein_complex_maps.protein_util as pu
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.plots.plot_bicluster as pbc
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.score_util as su
import protein_complex_maps.bicluster_generator as bg
import protein_complex_maps.annealer as anl

logging.basicConfig(level = logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s')


def main():

    parser = argparse.ArgumentParser(description="Non negative matrix foactorization for fractionation data")
    parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
                                            help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
                                            help="Protein ids in which to anaylze")
    parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
                                            help="Filename of output plot")
    parser.add_argument("--ignore_missing", action="store_true", dest="ignore_missing", required=False, default=False,
                                            help="Ignore missing protein ids in msds")
    parser.add_argument("--genenames", action="store_true", dest="genenames", required=False, default=False,
                                            help="Set labels to genenames")
    parser.add_argument("--rank", action="store", type=int, dest="rank", required=False, default=2,
                                            help="Rank or number of clusters")
    parser.add_argument("--max_iter", action="store", type=int, dest="max_iter", required=False, default=30,
                                            help="Maximum number of iterations for factorization")
    parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, default=0.0,
                                            help="Coefficient threshold for membership in a bicluster, default 0.0")

    args = parser.parse_args()

    msds = pickle.load( open( args.msds_filename, "rb" ) )

    if args.proteins != None:
            data_set, new_id_map = msds.get_subdata_matrix(args.proteins, ignoreNonExistingIds=args.ignore_missing) 
    else:
            data_set = msds.get_data_matrix()
            new_id_map = msds.get_name2index()


    if args.genenames:
            print "new_id_map.values: %s" % (new_id_map.values(),)
            genename_map = pu.get_genenames_uniprot( new_id_map.values() )
            print new_id_map
            print genename_map
            gene_id_map = dict()
            for i in xrange(len(data_set)):
                    print i
                    #print genename_map[new_id_map[i]]
                    #kdrew: sometimes no genename is returned for certain ids, default to original id
                    try:
                            if genename_map[new_id_map[i]] == None:
                                    gene_id_map[i] = new_id_map[i]
                            else:
                                    gene_id_map[i] = genename_map[new_id_map[i]]
                    except KeyError:
                            gene_id_map[i] = new_id_map[i]

            new_id_map = gene_id_map

    #logging.debug(clean_data_matrix)

    #kdrew: make all rows sum to 1.0
    #data_set = nu.normalize_over_columns(data_set)
    #kdrew: make all columns sum to 1.0
    data_set = nu.normalize_over_rows(data_set)

    V = data_set 

    print V

    # Run LSNMF algorithm
    # Returned object is fitted factorization model. Through it user can access quality and performance measures.
    # The fctr_res's attribute `fit` contains all the attributes of the factorization.  
    #fctr = nimfa.mf(V, method = "lsnmf", max_iter = args.max_iter, rank = args.rank)
    #fctr_res = nimfa.mf_run(fctr)

    fctr = nimfa.mf(V, method = "lsnmf", max_iter = args.max_iter, rank = args.rank )
    fctr_res = nimfa.mf_run(fctr)

    # Basis matrix.
    W = fctr_res.basis()
    print "Basis matrix"
    print W

    # Mixture matrix. 
    H = fctr_res.coef()
    print "Coef"
    #print H

    # Print the loss function according to Kullback-Leibler divergence. By default Euclidean metric is used.
    print "Distance Kullback-Leibler: %5.3e" % fctr_res.distance(metric = "kl")

    # Compute generic set of measures to evaluate the quality of the factorization
    sm = fctr_res.summary()
    # Print residual sum of squares (Hutchins, 2008). Can be used for estimating optimal factorization rank.
    print "Rss: %8.3f" % sm['rss']
    # Print explained variance.
    print "Evar: %8.3f" % sm['evar']
    # Print actual number of iterations performed
    print "Iterations: %d" % sm['n_iter']

    # Print estimate of target matrix V 
    #print "Estimate"
    #print np.dot(W, H)

    rsscore_obj = rsu.RandomSamplingScore(data_set, su.multiple_dot, sample_module=np.random)
    score_function = rsscore_obj.zscore_all
    bcgen = bg.BiclusterGenerator(rsscore_obj.zscore_all_neg, iterations=2500, starting_temperature = 0.00000001, random_module=np.random)
    ratequench_annealer = anl.RateQuenchAnnealer( bcgen.get_montecarlo(), quench_iteration=2400, scale=0.9, rate=0.3, recent_iterations=10, adjust_scale=True )
    bcgen.set_annealer(ratequench_annealer)


    bicluster1 = bc.Bicluster(rows=[0,1,2,3,4,5,14], cols=[55,56], random_module=random)
    print "lid zscore: %s" % (rsscore_obj.zscore_all(bicluster1.get_submatrix(data_set)),)
    print "lid cols zscore: %s" % (rsscore_obj.zscore_columns(data_set, bicluster1))
    print "lid rows zscore: %s" % (rsscore_obj.zscore_rows(data_set, bicluster1))

    print "Clusters"
    for i in range(args.rank):
        Wtransvec = np.transpose(W)[i]
        print "W",i, Wtransvec
        rows = np.ravel(np.where(Wtransvec > args.threshold)[1])
        Hvec = H[i]
        print "H",i, Hvec
        cols = np.ravel(np.where(Hvec > args.threshold)[1])
        print rows
        print cols
        bicluster1 = bc.Bicluster(rows=rows, cols=cols, random_module=random)
        pbc.plot_bicluster( data_set, bicluster1, savefilename="/home/kdrew/public_html/test/bicluster_test%s.pdf" % (i,), ylim_max=True)

        #kdrew: monte carlo optimize using nnmf cluster as seed
        bicluster_mc = bcgen.generator(data_set, bicluster1)
	eval_dict = bcgen.evaluate( data_set, len(bcgen.biclusters)-1 )

        try:
            print "zscore: %s" % (rsscore_obj.zscore_all(bicluster1.get_submatrix(data_set)),)
            print "cols zscore: %s" % (rsscore_obj.zscore_columns(data_set, bicluster1))
            print "rows zscore: %s" % (rsscore_obj.zscore_rows(data_set, bicluster1))
        except:
            continue

	for t in eval_dict.keys():
            logging.info("random %s, mean: %s, std: %s, zscore: %s" % ( t, eval_dict[t]['mean'], eval_dict[t]['std'], eval_dict[t]['zscore'] ))

	pbc.plot_bicluster(data_set, bicluster_mc, savefilename="/home/kdrew/public_html/test/bicluster_testMC%s.pdf" % (i), ylim_max=True)

    data_subplots = []
    f, data_subplots = plt.subplots(len(W),1,sharex='col')

    for i, data_r in enumerate(W):
        #print "pos: %s name: %s" % (i, new_id_map[i])
        data_row = np.array(data_r.reshape(-1))[0]
        barcolor = "gray"

        data_subplots[i].bar(np.arange(len(data_row)), data_row, align='center', facecolor=barcolor)
        data_subplots[i].axes.set_yticklabels([],visible=False)
        data_subplots[i].set_ylabel("%s" % (new_id_map[i],),rotation='horizontal', color=barcolor, fontsize=6, verticalalignment="center",horizontalalignment="right")


    plt.savefig(args.plot_filename)


if __name__ == "__main__":
	main()

