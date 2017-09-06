
import argparse
import numpy as np
from os.path import basename 


import protein_complex_maps.hhsuite.hhdomain as hhd

def main():

    parser = argparse.ArgumentParser(description="Generate domain list from hhresults by greedy selection")
    parser.add_argument("--complex_filename", action="store", dest="complex_filename", required=True, 
                                            help="Filenames of complexes")
    parser.add_argument("--hhr_filenames", action="store", dest="hhr_filenames", nargs='+', required=True, 
                                            help="Filenames of hhsuite results")
    parser.add_argument("--prob_threshold", action="store", type=float, dest="prob_threshold", required=False, default=90.0,
                                            help="Only consider domain hits above probability threshold (note: hhsuite prob is x100), default = 90.0")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                            help="Output of SCOP frequency table for each complex")
    parser.add_argument("--include_scop_levels", action="store_true", dest="include_scop_levels", required=False, default=False,
                                            help="Flag to include SCOP fold and superfamily levels in frequency counts")
    args = parser.parse_args()

    complexes = []
    complex_file = open(args.complex_filename,"rb")

    for line in complex_file.readlines():
        complexes.append(line.split())

    complex_file.close()

    domain_results = dict()

    #kdrew: generate nonoverlapping domain sets for each protein
    for hhr_filename in args.hhr_filenames:
        #print hhr_filename

        uniprot_id = basename(hhr_filename).split('_')[1]
        print uniprot_id

        greedy_list = hhd.greedy_domain_map(hhr_filename, args.prob_threshold)

        domain_results[uniprot_id] = greedy_list

        #for greedy_result in greedy_list:
        #    print "%s : %s-%s : %s" % (greedy_result['hit'], greedy_result['query_hmm_begin'], greedy_result['query_hmm_end'], greedy_result['prob'] )
        #
        #print ""


    #kdrew: loop through all complexes and make a list of every scop hit
    complexes_scop_dict = dict()
    for i, c in enumerate(complexes):
        print i
        complexes_scop_dict[i] = []
        for prot_id in c:
            try:
                print "%s: " % (prot_id,)
                for d in domain_results[prot_id]:
                    print "\t%s : %s-%s : %s" % (d['hit'], d['query_hmm_begin'], d['query_hmm_end'], d['prob'] )
                    #kdrew: parse scop id from hit record
                    scop_id = d['hit'].split()[1]
                    #kdrew: add scop id to complex list 
                    complexes_scop_dict[i].append(scop_id)

            except KeyError:
                print "%s: No results" % (prot_id,)

    #kdrew: create arrays of scop frequencies for each complex
    complexes_scop_freq_dict = dict()
    for i, c in enumerate(complexes):
        complexes_scop_freq_dict[i] = dict()
        for scop_id in complexes_scop_dict[i]:
            try:
                complexes_scop_freq_dict[i][scop_id] += 1
            except KeyError:
                complexes_scop_freq_dict[i][scop_id] = 1

            if args.include_scop_levels:
                #kdrew: add count for superfamily level
                sf_id = '.'.join(scop_id.split('.')[0:3])
                try:
                    complexes_scop_freq_dict[i][sf_id] += 1
                except KeyError:
                    complexes_scop_freq_dict[i][sf_id] = 1

                #kdrew: add count for fold level
                fold_id = '.'.join(scop_id.split('.')[0:2])
                try:
                    complexes_scop_freq_dict[i][fold_id] += 1
                except KeyError:
                    complexes_scop_freq_dict[i][fold_id] = 1


    #print complexes_scop_freq_dict
    fout = open(args.output_filename,"wb")
    for i in xrange(len(complexes)):
        output_str = "complex:%s\t" % (i,)
        for scop_id in complexes_scop_freq_dict[i]:
            output_str += "%s:%s\t" % (scop_id, complexes_scop_freq_dict[i][scop_id],)

        print output_str
        fout.write(output_str)
        fout.write("\n")

    fout.close()


if __name__ == "__main__":
	main()


