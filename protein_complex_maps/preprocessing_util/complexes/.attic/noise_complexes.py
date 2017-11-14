
import argparse
import numpy as np
import random
import itertools as it


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


def remove_subunits( complexes, totalNumberOfSubunits=1, totalFractionOfSubunits=0.1, fraction_flag=False, shuffle=True):
    #kdrew: make a copy of complexes and shuffle order of subunits
    shuffle_complexes = []
    for comp in complexes:
        c = list(comp[:])
        random.shuffle(c)
        shuffle_complexes.append(c)
    
    current_totalNumberOfSubunits = totalNumberOfSubunits
    if fraction_flag:
        #kdrew: calculate the total number of subunits to remove by finding sum of all subunits and multiply by total fraction parameter
        #kdrew: example, randomly remove 1/10 of subunits
        current_totalNumberOfSubunits = int(sum([len(c) for c in shuffle_complexes]) * totalFractionOfSubunits)


    #print "current_totalNumberOfSubunits: %s" % current_totalNumberOfSubunits
    i = 0
    while i < current_totalNumberOfSubunits: 
        complex_index = 0
        #kdrew: randomly pick complex and remove subunit
        if shuffle:
            complex_index = random.randint(0,len(shuffle_complexes)-1)

        comp = shuffle_complexes[complex_index]
        c1 = comp[1:]
        if len(c1) >= 2:
            #kdrew: replace complex for complex with removed subunit 
            shuffle_complexes[complex_index] = c1
            i += 1
        elif len(c1) == 1:
            #kdrew: if singleton do not add back to set of complexes and remove two subunits
            i += 2
            del shuffle_complexes[complex_index] 
        else:
            i += 1
            del shuffle_complexes[complex_index]


    return shuffle_complexes


#kdrew: randomly remove subunits from complexes
#kdrew: fraction is the fraction of complexes to remove subunits from 
#kdrew: numberOfSubunits is the number of subunits to remove from each complex
#kdrew: fractionOfSubunits is the fraction of subunits to remove from each complex, alternative to numberOfSubunits
#kdrew: fraction_flag is a boolean to turn on fractionOfSubunits
#kdrew: shuffle is a boolean that will shuffle the order of complexes before removing
#kdrew: default - remove one subunit from a random tenth of the complexes
def remove_subunits_by_complex( complexes, fraction=0.1, numberOfSubunits=1, fractionOfSubunits=0.1, fraction_flag=False, shuffle=True):
    #kdrew: make a copy of complexes and shuffle order
    shuffle_complexes = complexes[:]
    if shuffle:
        random.shuffle(shuffle_complexes)
    
    #kdrew: get index of complexes to remove subunits
    indexOfSplit = int(len(shuffle_complexes) * fraction)

    #kdrew: get list of complexes to remove subunits
    complexes2removesubunits = shuffle_complexes[:indexOfSplit]

    total_subunits_removed = 0
    removesubunit_complexes = []
    for comp in complexes2removesubunits:
        current_numberOfSubunits = numberOfSubunits
        if fraction_flag:
            current_numberOfSubunits = int(len(comp) * fractionOfSubunits)

        #kdrew: randomize subunits
        c = list(comp)[:]
        random.shuffle(c)
        c1 = c[current_numberOfSubunits:]

        total_subunits_removed += current_numberOfSubunits
        #kdrew: do not add singletons
        if len(c1) > 1:
            removesubunit_complexes.append(c1)
        else:
            #kdrew: add additional subunit to total_subunits_removed because we are removing the singleton as well
            total_subunits_removed += 1

    final_complexes = removesubunit_complexes + shuffle_complexes[indexOfSplit:]

    return final_complexes, total_subunits_removed



def add_random_complexes( complexes, fraction=0.1 ):
    ids = set()
    for c in complexes:
        ids = ids.union(set(c)) 

    random_complexes = []
    ids_list = list(ids)
    numberOfComplexes = int(len(complexes)*fraction)
    for i in range(numberOfComplexes):
        random.shuffle(ids_list)
        #kdrew: grab random complex and use that size as the size of the new random complex
        sizeOfComplex = len(complexes[random.randint(0,len(complexes)-1)])
        random_complexes.append(ids_list[:sizeOfComplex])


    return complexes + random_complexes



#kdrew: recipe from https://docs.python.org/2/library/itertools.html
def powerset(iterable):
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) for r in range(len(s)+1))

def powerset_complexes( complexes ):
    ids = set()
    for c in complexes:
        ids = ids.union(set(c)) 
    full_powerset = powerset(ids)

    #kdrew: remove singletons
    full_powerset_filtered_list = [c for c in full_powerset if len(c) >= 2]
    return full_powerset_filtered_list

#kdrew: order complexes by size and remove complexes in large-size portion (or small-size portion depending remove_upper flag)
def remove_fraction_complexes( complexes, fraction = 0.5, remove_upper=True ):
    #kdrew: order complexes by size
    sorted_complexes = sorted(complexes, key=len)
    index = int( len(sorted_complexes) * fraction)
    if remove_upper:
        return sorted_complexes[:index]
    else:
        return sorted_complexes[index:]



LARGE_INT = 1000000
def remove_by_size_complexes( complexes, upper_size_threshold = LARGE_INT, lower_size_threshold = 1):
    filtered_complexes = [ c for c in complexes if len(c) > lower_size_threshold and len(c) <= upper_size_threshold ]
    return filtered_complexes

def breakup_complexes(complexes, breakup_fraction=0.1, keep_original=False ):
    
    #kdrew: make a copy of complexes and shuffle order
    shuffle_complexes = complexes[:]
    random.shuffle(shuffle_complexes)
    
    #kdrew: get index of complexes to breakup
    indexOfBreakupIds = int(len(shuffle_complexes) * breakup_fraction)

    #kdrew: get list of complexes to breakup
    complexes2breakup = shuffle_complexes[:indexOfBreakupIds]

    broken_complexes = []
    for comp in complexes2breakup:
        #kdrew: copy complex
        c = list(comp)[:]

        #kdrew: shuffle complex members
        random.shuffle(c)

        #kdrew: get index of half way point in complex list
        split_id = len(c)/2

        #kdrew: split complex
        c1 = c[:split_id]
        c2 = c[split_id:]

        #kdrew: do not add singletons
        if len(c1) > 1:
            broken_complexes.append(c1)
        if len(c2) > 1:
            broken_complexes.append(c2)

    final_complexes = broken_complexes + shuffle_complexes[indexOfBreakupIds:]

    if keep_original:
        final_complexes = final_complexes + shuffle_complexes[:indexOfBreakupIds]

    return final_complexes


def shuffle_complexes(complexes, shuffle_fraction=0.1 ):

    complex_id_list = []
    protein_id_list = []
    for i, cmplex in enumerate(complexes):
        #print i
        #print cmplex
        for prot in cmplex:
            complex_id_list.append(i)
            protein_id_list.append(prot)

    index_list = range(len(complex_id_list))

    #print complex_id_list
    #print protein_id_list
    #print index_list

    #kdrew: shuffle to get a random list of protein complex membership indices
    random.shuffle(index_list)
    #print index_list

    #kdrew: index of the shuffled index_list of shuffle_fraction
    indexOfShuffledIds = int(len(index_list) * shuffle_fraction)
    #print indexOfShuffledIds

    #kdrew: get top ids
    indices2shuffle = index_list[:indexOfShuffledIds]
    #print indices2shuffle

    #kdrew: shuffle top ids complex membership
    top_complex_id_list = [complex_id_list[x] for x in indices2shuffle]
    #print top_complex_id_list
    random.shuffle(top_complex_id_list)
    #print top_complex_id_list

    #kdrew: update the complex membership with the shuffled memberships
    for i,x in enumerate(indices2shuffle):
        complex_id_list[x] = top_complex_id_list[i]

    #print protein_id_list
    #print complex_id_list

    shuffled_complexes_dict = dict()
    for i, prot_id in enumerate(protein_id_list):
        try:
            shuffled_complexes_dict[complex_id_list[i]].add(prot_id)
        except KeyError:
            shuffled_complexes_dict[complex_id_list[i]] = set([prot_id])


    #kdrew: return list of complexes which are lists themselves
    shuffled_complexes = [list(xx) for xx in shuffled_complexes_dict.values()]
    return shuffled_complexes


if __name__ == "__main__":
	main()

