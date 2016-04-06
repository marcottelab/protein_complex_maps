
import argparse
import numpy as np

import protein_complex_maps.hhsuite.hhresults_parser as hrp

def main():

    parser = argparse.ArgumentParser(description="Generate domain list from hhresults by greedy selection")
    parser.add_argument("--hhr_filenames", action="store", dest="hhr_filenames", nargs='+', required=True, 
                                            help="Filenames of hhsuite results")
    parser.add_argument("--prob_threshold", action="store", type=float, dest="prob_threshold", required=False, default=90.0,
                                            help="Only consider domain hits above probability threshold (note: hhsuite prob is x100), default = 90.0")
    args = parser.parse_args()


    for hhr_filename in args.hhr_filenames:

        print hhr_filename

        greedy_list = greedy_domain_map(hhr_filename, args.prob_threshold)

        for greedy_result in greedy_list:
            print "%s : %s-%s : %s" % (greedy_result['hit'], greedy_result['query_hmm_begin'], greedy_result['query_hmm_end'], greedy_result['prob'] )

        print ""


def greedy_domain_map(hhr_filename, prob_threshold):

    results_list = hrp.parse_hhr_hitlist(hhr_filename)
    greedy_list = []

    #kdrew: loop through all results (assume sorted by probability)
    for result in results_list:
        #kdrew: check to see if result is above probablity threshold, if not break
        if result['prob'] < prob_threshold:
            break

        nooverlap = True
        #kdrew: compare result with all results already in the domain set
        for greedy_result in greedy_list:
            #kdrew: if results overlap flip nooverlap flag
            if hits_overlap(greedy_result, result):
                nooverlap = False
                break
        #kdrew: if there is no overlap between result and all other no overlapping results, store
        if nooverlap:
            greedy_list.append(result)

    return greedy_list



def hits_overlap(result1, result2):
    if result1['query_hmm_begin'] <= result2['query_hmm_end'] and result1['query_hmm_end'] >= result2['query_hmm_begin']:
        return True
    else:
        return False


if __name__ == "__main__":
	main()


