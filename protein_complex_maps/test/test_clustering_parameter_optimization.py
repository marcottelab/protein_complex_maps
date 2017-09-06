#!/usr/bin/python

# Unit tests for new clustering analysis functionality

import unittest
import protein_complex_maps.features.clustering_parameter_optimization as cpo

import numpy as np

class Object(object):
    pass


class ClusterTest(unittest.TestCase):

    def setUp(self,):

        input_network_list = ['A\tB\t1.0\n','B\tC\t1.0\n','A\tC\t1.0\n','A\tD\t0.2\n','C\tD\t0.2\n','D\tE\t0.7\n','E\tF\t0.8\n','D\tF\t0.75\n'] 
        ppi_scores = dict()
        for ppi in input_network_list:
            ppi_scores[frozenset([ppi.split()[0],ppi.split()[1]])] = float(ppi.split()[2])

        fraction = 1.0

        sizeOfTopNetwork = int(len(input_network_list)*float(fraction))
        network_list = input_network_list[0:sizeOfTopNetwork]
        threshold_score = float(network_list[-1].split()[2])
    

        args = Object()
        args.cfinder_license = '/home/kdrew/programs/CFinder-2.0.6--1448/licence.txt'
        args.cfinder_exe = '/home/kdrew/programs/CFinder-2.0.6--1448/CFinder_commandline64'
        args.temp_dir = '/tmp/'
        args.mcl_bin = 'mcl'
        args.clustone_jar = '/home/kdrew/programs/clusterone/cluster_one-1.0.jar'

        self.parameter_dict = dict()
        self.parameter_dict['network_list'] = network_list
        self.parameter_dict['ppi_scores'] = ppi_scores
        self.parameter_dict['args'] = args
        self.parameter_dict['size'] = str(2)
        self.parameter_dict['density'] = str(0.4)
        self.parameter_dict['overlap'] = str(0.8)
        self.parameter_dict['seed_method'] = 'nodes'
        self.parameter_dict['fraction'] = str(fraction)
        self.parameter_dict['threshold_score'] = threshold_score
        self.parameter_dict['i'] = '00'
        self.parameter_dict['inflation'] = '2.0'
        self.parameter_dict['cliquesize'] = '3'
        self.parameter_dict['timeout'] = '10'
        self.parameter_dict['twostep_combination'] = ['clusterone','mcl']

    def testClustering(self,):

        predicted_clusters, i  = cpo.cluster_helper(self.parameter_dict)
        print predicted_clusters
        #assert (Y[0][0] == 4)
        #assert (Y[0][1] == 7)
        #assert (Y[1][0] == 8)
        #assert (Y[1][1] == 9)

    def testClustering2(self,):

        self.parameter_dict['twostep_combination'] = ['cfinder','agglomod']
        predicted_clusters, i  = cpo.cluster_helper(self.parameter_dict)
        print predicted_clusters


if __name__ == "__main__":
    unittest.main()


