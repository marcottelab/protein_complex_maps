
import argparse
import numpy as np
import random

import protein_complex_maps.clustering.clustering_parameter_optimization as cpo


def main():

    parser = argparse.ArgumentParser(description="Trim subunits from complexes that do not meet edge weight threshold")
    parser.add_argument("--input_complexes", action="store", dest="input_complexes", required=True, 
                                            help="Filename of input complexes (one complex per line, ids space/tab separated)")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True, 
                                            help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--threshold", action="store", dest="threshold", type=float, required=True, 
                                            help="Trim complex subunits that do not have any edges with weight greater than threshold")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Filename of where the output should go")
    args = parser.parse_args()

    #kdrew: read in the input network into a string
    with open (args.input_network, "r") as input_network_file:
        input_network_list = input_network_file.readlines()

    ppi_scores = dict()
    for ppi in input_network_list:
        ppi_scores[frozenset([ppi.split()[0],ppi.split()[1]])] = float(ppi.split()[2])

    complexes = []
    f = open(args.input_complexes,"rb")
    for line in f.readlines():
        complexes.append(line.split())
    f.close()

    trimed_complexes = cpo.trim_clusters2threshold(complexes, args.threshold, ppi_scores)

    outfile = open(args.output_filename,"wb")
    for c in trimed_complexes:
        outfile.write("%s\n" % (" ".join(c)))

    outfile.close()

    
if __name__ == "__main__":
	main()

