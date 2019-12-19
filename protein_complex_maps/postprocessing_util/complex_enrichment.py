    
import numpy as np
import argparse
import itertools as it

import requests
import subprocess as sp

import pandas as pd

try:
   import cPickle as pickle
except:
   import pickle


def main():

    parser = argparse.ArgumentParser(description="Calculate term enrichment for complexes using gprofiler")
    parser.add_argument("--complexes", action="store", dest="complex_filename", required=True, 
                                            help="Filename of complexes, format one cluster per line, ids space separated")
    parser.add_argument("--background", action="store", dest="background_filename", required=False, default=None,
                                            help="Filename of background ids, format one id per line")
    #parser.add_argument("--secondary_complexes", action="store", dest="secondary_complex_filename", required=True, 
    #                                        help="Filename of secondary complexes, format one cluster per line, ids space separated")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=False, default=None, 
                                            help="Filename to output difference complexes")
    parser.add_argument("--correction_method", action="store", dest="correction_method", required=False, default='fdr', 
                                            help="Correction method for multiple hypothesis testing {gSCS,fdr,bonferroni}, default = fdr")
    parser.add_argument("--resume_from_checkpoint", action="store_true", dest="resume_from_checkpoint", required=False, default=False, 
                                            help="Option to resume from checkpointed file: results_checkpoint.pkl, default = False")

    args = parser.parse_args()

    #ids = "Q13772 Q13200 O00487 Q9BRP4 O43242 P62195 P43686 Q99460 P62333 P35998 P62191 Q16186 P17980 Q05086 Q15008 Q9Y5K5 Q9UNM6 O00231 O00232 Q53HC0 O75832 P51665 P48556 P55036 Q16401 O00233 Q96EN9"

    #if args.output_filename != None:
    #    output_file = open(args.output_filename,'wb')

    background_proteins = []
    if args.background_filename != None:
        background_file = open(args.background_filename,'rb')
        for line in background_file.readlines():
            background_proteins.append(line.split()[0])

    complex_file = open(args.complex_filename,'rb')
    results_df = None
    index_count = 0

    if args.resume_from_checkpoint:
            results_pickle = open('results_checkpoint.pkl', 'rb')
            results_df = pickle.load(results_pickle)
            results_pickle.close()

            index_count = max(results_df['index']) + 1
            checkpoint_complex_id = max(results_df['complex_id'])

    for i, complex_line in enumerate(complex_file.readlines()):

        print "complex: %s" % (i,)
        #kdrew: skip complexes if already from checkpoint
        if args.resume_from_checkpoint:
            if i <= checkpoint_complex_id:
                print "skipping: resuming from checkpoint file"
                continue
        
        if len(complex_line.split()) == 0:
            continue

        complex_query_ids = list(set(complex_line.split()))
        r = requests.post(
                url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                json={
                    'organism':'hsapiens',
                    'query':complex_query_ids,
                    'no_iea':True,
                    'no_evidences':False,
                    'domain_scope':'custom',
                    'background':background_proteins,
                    }
                )
        #kdrew: does not return list in same order as input
        query_ids = r.json()['meta']['genes_metadata']['query']['query_1']['mapping'].keys()
        print query_ids
        ordered_query_ids = [x for x in complex_query_ids if x in query_ids]
        for res in r.json()['result']:
            res['index'] = index_count
            index_count += 1
            res['complex_id'] = i
            print res
            #kdrew: this field makes converting to a dataframe difficult
            #res['intersections'] = ' '.join(res['intersections'])
            #kdrew: convert so that we know which genes are annotated
            res['intersection_genes'] = ' '.join([yy for j,yy in enumerate(ordered_query_ids) if len(res['intersections'][j]) > 0])
            #kdrew: combine parents into a single string otherwise it gets separated into two entries, 
            #kdrew: put string back into list because pandas complains about all scalars and no index but can't set index through from_dict function
            res['parents'] = [' '.join(res['parents'])]
            del res['intersections']
            #kdrew: pandas to_csv chokes on writing the description of an enrichment due to not unicode encoding
            #kdrew: UnicodeEncodeError: 'ascii' codec can't encode character u'\xef' in position 572: ordinal not in range(128)
            res['description'] = res['description'].encode('utf-8')
            if results_df is None:
                results_df = pd.DataFrame.from_dict(res)
            else:
                results_df = pd.concat([results_df, pd.DataFrame.from_dict(res)])

        #kdrew: checkpoint just in case of failures
        if i % 500 == 0:
            results_pickle = open('results_checkpoint.pkl', 'wb')
            pickle.dump(results_df, results_pickle)
            results_pickle.close()

    results_df = results_df.set_index("index")
    print results_df
    #kdrew: pickle just in case pandas fails again
    results_pickle = open('results.pkl', 'wb')
    pickle.dump(results_df, results_pickle)
    
    results_df.to_csv(args.output_filename, encoding='utf-8')


if __name__ == "__main__":
        main()

