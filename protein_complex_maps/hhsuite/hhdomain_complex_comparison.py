
import argparse
import numpy as np
import operator
import pickle as p

import protein_complex_maps.hhsuite.hhdomain as hhd

INF_ZSCORE = 100.0

def main():

    parser = argparse.ArgumentParser(description="Compare complexes based on SCOP frequencies")
    parser.add_argument("--hhdomain_scop_filename", action="store", dest="hhdomain_scop_filename", required=True, 
                                            help="Filename of lists of scop domains per protein")
    parser.add_argument("--complex_filename", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes")
    parser.add_argument("--shuffle_complex_filenames", nargs='+',  action="store", dest="shuffle_filenames", required=True, 
                                            help="Filenames of shuffled complex scop frequency tables")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                            help="Filename to store scoring output between clusters")
    parser.add_argument("--protein_edge_output_filename", action="store", dest="protein_edge_output_filename", required=True, 
                                            help="Filename to store output protein edges")
    parser.add_argument("--results_pickle_filename", action="store", dest="results_pickle_filename", required=False, default=None,
                                            help="Filename to store results in pickle file")
    args = parser.parse_args()

    complexes = []
    complex_file = open(args.complex_filename,"rb")
    for line in complex_file.readlines():
        complexes.append(line.split())
    complex_file.close()

    hhdomain_dict = dict()
    hhdomain_file = open(args.hhdomain_scop_filename,"rb")
    for line in hhdomain_file.readlines():
        protid = line.split(',')[0]
        hhdomain_dict[protid] = [x.rstrip() for x in line.split(',')[1:]]
    hhdomain_file.close()

    shuffled_complexes = dict()
    for shuffle_filename in args.shuffle_filenames:
        shuffle_complexes = []
        shuffle_file = open(shuffle_filename, "rb")
        for line in shuffle_file.readlines():
            shuffle_complexes.append(line.split())
        shuffle_file.close()
        shuffled_complexes[shuffle_filename] = shuffle_complexes

    comparison_results = dict()
    #kdrew: compare every complex to every other complex
    for i, c1 in enumerate(complexes):
        for j, c2 in enumerate(complexes):
            #kdrew: overlap_score is the number of common scop ids, overlap_domains is a list of ((protein_id, scop_id, iteration), (protein_id2, scop_id, iteration)) tuples
            overlap_score, overlap_domains, protein_count = complex_comparison(c1, c2, hhdomain_dict)
            #kdrew: if there is any overlap
            if overlap_score > 0:
                #kdrew: compare first complex to random set of complexes of the same size as second complex, store overlap_scores in index 0
                rand_scores = [complex_comparison(c1,rand1[j], hhdomain_dict)[0] for rand1 in shuffled_complexes.values()]
                zscore = (overlap_score - np.mean(rand_scores)) / np.std(rand_scores)
                #kdrew: some of the zscores evaluate to 'inf', setting to a large number ex 100, didn't want to make it too large because the zscores are the widths in cytoscape
                if np.isinf(zscore):
                    zscore = INF_ZSCORE
                comparison_results[(i,j)] = {'overlap_score':overlap_score, 'overlap_domains':overlap_domains, 'zscore':zscore, 'protein_count':protein_count}

    if args.results_pickle_filename != None:
        p.dump(comparison_results,open(args.results_pickle_filename,"wb"))

    fout = open(args.output_filename,"wb")
    prot_out = open(args.protein_edge_output_filename,"wb")
    fout.write("i, j, overlap_score, zscore, jaccard, protein_count\n") 
    for i,j in comparison_results.keys():
        if comparison_results[(i,j)]['overlap_score'] > 2:
        
            fout.write("%s,%s,%s,%s,%s,%s\n" % (i,j,comparison_results[(i,j)]['overlap_score'],
                                                            comparison_results[(i,j)]['zscore'], 
                                                            jaccard_index(complexes[i],complexes[j]), 
                                                            comparison_results[(i,j)]['protein_count'],))
            #print "%s , %s, overlap_score: %s, zscore: %s, jaccard: %s, protein_count: %s" % (i,j,comparison_results[(i,j)]['overlap_score'],
                                                                                                    #comparison_results[(i,j)]['zscore'], 
                                                                                                    #jaccard_index(complexes[i],complexes[j]), 
                                                                                                    #comparison_results[(i,j)]['protein_count'],)

            for domain_pair in comparison_results[(i,j)]['overlap_domains']:
                #print "%s (%s) %s" % (domain_pair[0][0], domain_pair[0][1], domain_pair[1][0], )
                #kdrew: output is in cytoscape edge format
                prot_out.write("%s_%s (%s_%s) %s_%s\n" % (i, domain_pair[0][0], domain_pair[0][1], domain_pair[0][2], j, domain_pair[1][0], ))

    fout.close()
    prot_out.close()


def complex_comparison(complex1, complex2, hhdomain_dict):
    c1_scop_set = complex_scop_set(complex1,hhdomain_dict)
    #print c1_scop_set
    c2_scop_set = complex_scop_set(complex2,hhdomain_dict)
    #print c2_scop_set
    
    overlap_domains = set()
    protein_set = set()
    for domain in c1_scop_set:
        #kdrew: copy scop set so we can remove items and not have issues with iterating while changing the set
        c2_scop_set_copy = set([x for x in c2_scop_set])
        for domain2 in c2_scop_set_copy:
            #kdrew: if scop_ids match and iteration match
            if domain[1] == domain2[1] and domain[2] == domain2[2]:
                overlap_domains.add((domain,domain2))
                #kdrew: add proteins in match to set
                protein_set.add(domain[0])
                protein_set.add(domain2[0])
                #kdrew: remove domain2 from set so we don't try and match it again
                c2_scop_set.remove(domain2)
                break

    return len(overlap_domains), overlap_domains, len(protein_set)

def complex_scop_set(c, hhdomain_dict):
    scop_set = set()
    for protid in c:
        try:
            domains = hhdomain_dict[protid]
            for scop_i in set(domains):
                if scop_i != '':
                    for i in range(hhdomain_dict[protid].count(scop_i)):
                        scop_set.add((protid,scop_i,i))
        except KeyError:
            continue
    return scop_set

def protein_comparison(prot1, prot2, hhdomain_dict):
    scop_intersection = dict()
    for scop_i in hhdomain_dict[prot1]:
        if scop_i != '' and scop_i in hhdomain_dict[prot2]:
            scop_intersection[scop_i] = min(hhdomain_dict[prot1].count(scop_i), hhdomain_dict[prot2].count(scop_i))

    #kdrew: I think this is unnecessary, if it didn't match in prot1 there won't be a match in prot2
    #for scop_j in hhdomain_dict[prot2]:
    #    if scop_j != '' and scop_j not in scop_intersection.keys() and scop_j in hhdomain_dict[prot1]:
    #        scop_intersection[scop_j] = min(hhdomain_dict[prot1].count(scop_j), hhdomain_dict[prot2].count(scop_j))

    return scop_intersection



def temp_func():
    complexes_scop_freq_dict = dict()
    freq_file = open(args.scop_freq_filename,"rb")
    for line in freq_file.readlines():
        i = line.split()[0].split(':')[1]
        complexes_scop_freq_dict[i] = dict()
        for sid_freq in line.split()[1:]:
            sid = sid_freq.split(':')[0]
            freq = int(sid_freq.split(':')[1])
            complexes_scop_freq_dict[i][sid] = freq

    freq_file.close()

    shuffled_scop_freq_list = []
    for shuffle_file in args.shuffle_scop_freq_filenames:
        shuffled_complexes_scop_freq_dict = dict()
        freq_file = open(shuffle_file,"rb")
        for line in freq_file.readlines():
            i = line.split()[0].split(':')[1]
            shuffled_complexes_scop_freq_dict[i] = dict()
            for sid_freq in line.split()[1:]:
                sid = sid_freq.split(':')[0]
                freq = int(sid_freq.split(':')[1])
                shuffled_complexes_scop_freq_dict[i][sid] = freq

        freq_file.close()
        shuffled_scop_freq_list.append(shuffled_complexes_scop_freq_dict)




    complex_similarity = dict()
    complex_similarity_zscore = dict()
    #kdrew: for every complex vs every other complex, compute similarity score (dot product or correlation) of scop frequency tables
    for i in complexes_scop_freq_dict:
        complex_similarity[i] = dict()
        complex_similarity_zscore[i] = dict()
        for j in complexes_scop_freq_dict:
            similarity_metric = min_sum_score(complexes_scop_freq_dict[i], complexes_scop_freq_dict[j])
            #print "%s,%s : %s" % (i,j,similarity_metric)
            if similarity_metric > 0.0:
                shuffled_similarity_list = []
                for ii in range(len(shuffled_scop_freq_list)):
                    shuffled_similarity_list.append(min_sum_score(complexes_scop_freq_dict[i], shuffled_scop_freq_list[ii][j]))


                complex_similarity[i][j] = similarity_metric
                #kdrew: TODO: need to check for divide by zero
                complex_similarity_zscore[i][j] = (1.0*similarity_metric - np.mean(shuffled_similarity_list)) / np.std(shuffled_similarity_list)

    for i in complex_similarity_zscore:
        sorted_x = sorted(complex_similarity_zscore[i].items(), key=operator.itemgetter(1))
        #print "complex: %s" % (i,)
        for j,zscore in sorted_x:
            print "%s,%s,%s,%s,%s" % (i, j, zscore, complex_similarity[i][j], jaccard_index(complexes[int(i)],complexes[int(j)]), )

def jaccard_index(x, y):
    sx = set(x)
    sy = set(y)
    return 1.0 * len(sx.intersection(sy)) / len(sx.union(sy))


def min_sum_score(complex_freq1, complex_freq2):
    similarity_metric = 0.0
    for scop_i in complex_freq1:
        try:
            similarity_metric += min(complex_freq2[scop_i],complex_freq1[scop_i])
        except KeyError:
            continue
    for scop_j in complex_freq2:
        if scop_j not in complex_freq1.keys():
            try:
                similarity_metric += min(complex_freq2[scop_j],complex_freq1[scop_j])
            except KeyError:
                continue

    return similarity_metric
        
def dot_product(complex_freq1, complex_freq2):
    similarity_metric = 0.0
    for scop_i in complex_freq1:
        try:
            similarity_metric += complex_freq2[scop_i] * complex_freq1[scop_i]
        except KeyError:
            continue
    for scop_j in complex_freq2:
        if scop_j not in complex_freq1.keys():
            try:
                similarity_metric += complex_freq2[scop_j] * complex_freq1[scop_j]
            except KeyError:
                continue

    return similarity_metric

if __name__ == "__main__":
	main()


