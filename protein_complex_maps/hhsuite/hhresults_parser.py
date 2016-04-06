
import argparse
import numpy as np


def main():

    parser = argparse.ArgumentParser(description="Parses hit list from hhsuite result files (hhr format)")
    parser.add_argument("--hhr_filenames", action="store", dest="hhr_filenames", nargs='+', required=True, 
                                            help="Filenames of hhsuite results")
    args = parser.parse_args()


    hhr_dict = dict()
    for hhr_filename in args.hhr_filenames:
        hhr_dict[hhr_filename] = parse_hhr_hitlist(hhr_filename)
        print hhr_dict[hhr_filename]


def parse_hhr_hitlist(hhr_filename):
    f = open(hhr_filename,"rb")

    in_results_summary = False
    results_list = []

    for line in f.readlines():
        if "No Hit" in line :
            in_results_summary = True
            continue

        if in_results_summary and len(line) > 1:
            in_results_summary = True
            #print line
            hhresult_dict = dict()
            hhresult_dict['no'] = int(line[0:3])
            hhresult_dict['hit'] = line[4:34]
            hhresult_dict['prob'] = float(line[35:40])
            hhresult_dict['evalue'] = float(line[41:48])
            hhresult_dict['pvalue'] = float(line[49:56])
            hhresult_dict['score'] = float(line[58:65])
            hhresult_dict['ss'] = float(line[66:70])
            hhresult_dict['cols'] = int(line[70:75])
            hhresult_dict['query_hmm'] = line[75:85]
            hhresult_dict['query_hmm_begin'] = int(hhresult_dict['query_hmm'].split('-')[0])
            hhresult_dict['query_hmm_end'] = int(hhresult_dict['query_hmm'].split('-')[1])
            hhresult_dict['template_hmm'] = line[85:100]
            hhresult_dict['template_hmm_begin'] = int(hhresult_dict['template_hmm'].split('(')[0].split('-')[0])
            hhresult_dict['template_hmm_end'] = int(hhresult_dict['template_hmm'].split('(')[0].split('-')[1])
            hhresult_dict['template_hmm_match_states'] = int(hhresult_dict['template_hmm'].split('(')[1].translate(None, '()'))

            #print hhresult_dict
            results_list.append(hhresult_dict)

        else:
            in_results_summary = False

    f.close()
    return results_list    



if __name__ == "__main__":
	main()


