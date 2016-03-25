
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

    args = parser.parse_args()


    graph = nx.read_edgelist(args.input_network, nodetype=str, data=(('weight',float),))

    newman = ag.NewmanGreedy(graph)
    print newman.quality_history
    print newman.get_clusters()




if __name__ == "__main__":
    main()


