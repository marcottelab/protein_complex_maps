
import argparse
import numpy as np
import random


def main():

    parser = argparse.ArgumentParser(description="Add noise to complexes (add and remove complex membership)")
    parser.add_argument("--input_complexes", action="store", dest="input_complexes", required=True, 
                                            help="Filename of input complexes (one complex per line, ids space/tab separated)")
    parser.add_argument("--shuffle_fraction", action="store", dest="shuffle_fraction", type=float, required=True, 
                                            help="Fraction of complex memberships to shuffle, range: 0.0 - 1.0")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                                            help="Filename of where the output should go")
    args = parser.parse_args()



    complexes = []
    f = open(args.input_complexes,"rb")
    for line in f.readlines():
        complexes.append(line.split())
    f.close()

    noise_complexes = shuffle_complexes(complexes, args.shuffle_fraction)
    
    outfile = open(args.output_filename,"wb")
    for c in noise_complexes:
        outfile.write("%s\n" % (" ".join(c)))

    outfile.close()


def shuffle_complexes(complexes, shuffle_fraction=0.1 ):

    complex_id_list = []
    protein_id_list = []
    for i, cmplex in enumerate(complexes):
        for prot in cmplex:
            complex_id_list.append(i)
            protein_id_list.append(prot)

    index_list = range(len(complex_id_list))

    print complex_id_list
    print protein_id_list
    print index_list

    #kdrew: shuffle to get a random list of protein complex membership indices
    random.shuffle(index_list)
    print index_list

    #kdrew: index of the shuffled index_list of shuffle_fraction
    indexOfShuffledIds = int(len(index_list) * shuffle_fraction)
    print indexOfShuffledIds

    #kdrew: get top ids
    indices2shuffle = index_list[:indexOfShuffledIds]
    print indices2shuffle

    #kdrew: shuffle top ids complex membership
    top_complex_id_list = [complex_id_list[x] for x in indices2shuffle]
    print top_complex_id_list
    random.shuffle(top_complex_id_list)
    print top_complex_id_list

    #kdrew: update the complex membership with the shuffled memberships
    for i,x in enumerate(indices2shuffle):
        complex_id_list[x] = top_complex_id_list[i]

    print protein_id_list
    print complex_id_list

    shuffled_complexes_dict = dict()
    for i, prot_id in enumerate(protein_id_list):
        try:
            shuffled_complexes_dict[complex_id_list[i]].add(prot_id)
        except KeyError:
            shuffled_complexes_dict[complex_id_list[i]] = set([prot_id])


    #kdrew: return list of complexes which are lists themselves
    shuffled_complexes = [list(xx) for xx in shuffled_complexes_dict.values()]
    return shuffled_complexes


    #complexes_dict = dict()
    #
    #for i, cmplex in enumerate(complexes):
    #    for prot in cmplex:
    #        try:
    #            complexes_dict[prot].append(i)
    #        except KeyError:
    #            complexes_dict[prot] = [i]
    #


    


if __name__ == "__main__":
	main()

