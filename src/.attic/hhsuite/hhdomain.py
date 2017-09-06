
import argparse
import numpy as np

from os.path import basename

import protein_complex_maps.hhsuite.hhresults_parser as hrp
import protein_complex_maps.protein_util as pu

def main():

    parser = argparse.ArgumentParser(description="Generate domain list from hhresults by greedy selection")
    parser.add_argument("--hhr_filenames", action="store", dest="hhr_filenames", nargs='+', required=True, 
                                            help="Filenames of hhsuite results")
    parser.add_argument("--prob_threshold", action="store", type=float, dest="prob_threshold", required=False, default=90.0,
                                            help="Only consider domain hits above probability threshold (note: hhsuite prob is x100), default = 90.0")
    parser.add_argument("--parse_uniprot_from_filename", action="store_true", dest="parse_uniprot_from_filename", required=False, default=False,
                                            help="Parse uniprot id from hhr filename")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                            help="Filename to write output to")

    parser.add_argument("--convert_id", action="store_true", dest="convert_id", required=False, default=False,
                                            help="Convert id from one identifier to another")
    parser.add_argument("--to_id", action="store", dest="to_id", required=False, default="ACC", 
                                            help="convert ids to this type, (P_ENTREZGENEID,ACC,etc), default=ACC")
    parser.add_argument("--reviewed", action="store_true", dest="reviewed", required=False, default=False,
                                            help="map only to reviewed ids, default=False")
    args = parser.parse_args()

    

    domain_list_dict = dict()
    for hhr_filename in args.hhr_filenames:
        greedy_list = greedy_domain_map(hhr_filename, args.prob_threshold)
        domain_list_dict[hhr_filename] = [x['hit'].split()[1] for x in greedy_list]

    uniprot_domain_list_dict = dict()
    if args.parse_uniprot_from_filename:
        for hhr_filename in domain_list_dict:
            uniprot_id = basename(hhr_filename).split('_')[1]
            uniprot_domain_list_dict[uniprot_id] = domain_list_dict[hhr_filename]
        domain_list_dict = uniprot_domain_list_dict

        if args.convert_id:
            ACC2outputID_map = pu.map_protein_ids(uniprot_domain_list_dict.keys(), "ACC", args.to_id, reviewed=args.reviewed)
            convert_id_domain_list_dict = dict()
            for uniprot_id in uniprot_domain_list_dict:
                try:
                    converted_id = ACC2outputID_map[uniprot_id][0]
                    convert_id_domain_list_dict[converted_id] = uniprot_domain_list_dict[uniprot_id]
                except KeyError:
                    print "cannot find uniprot_id: %s" % (uniprot_id)
                    convert_id_domain_list_dict[uniprot_id] = uniprot_domain_list_dict[uniprot_id]

            domain_list_dict = convert_id_domain_list_dict

    fout = open(args.output_filename,"wb")
    for protid in domain_list_dict:
        out_str = "%s,%s\n" % (protid,','.join(domain_list_dict[protid]))
        fout.write(out_str)

    fout.close()



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


