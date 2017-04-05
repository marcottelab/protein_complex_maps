
import os
import argparse
import pickle
import random
import numpy as np
import pandas as pd
import itertools as it
import subprocess as sp
import multiprocessing as mp
import tempfile as tf

import networkx as nx
import agglomcluster.agglomod as ag

def main():

    parser = argparse.ArgumentParser(description="Runs Agglomerative clustering")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True, 
                                    help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                    help="Filename to output clusters")

    args = parser.parse_args()


    graph = nx.read_edgelist(args.input_network, nodetype=str, data=(('weight',float),))

    newman = ag.NewmanGreedy(graph)
    #print newman.quality_history
    fout = open(args.output_filename,"wb")
    for clust in newman.get_clusters():
        fout.write("%s\n" % " ".join(clust))

    fout.close()




if __name__ == "__main__":
    main()


