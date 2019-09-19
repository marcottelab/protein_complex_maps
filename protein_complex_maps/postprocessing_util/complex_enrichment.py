    
import numpy as np
import argparse
import itertools as it

import requests
import subprocess as sp

import pandas as pd


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
    for i, complex_line in enumerate(complex_file.readlines()):
        #if args.output_filename != None:
        #    output_file.write("complex: %s\n" % (i,))
        #    output_file.write("""#  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list\n""")

        print "complex: %s" % (i,)
        #print """#  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list"""
        

        r = requests.post(
                url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
                json={
                    'organism':'hsapiens',
                    'query':complex_line.split(),
                    'no_iea':True,
                    'no_evidences':False,
                    'domain_scope':'custom',
                    'background':background_proteins,
                    }
                )
        for res in r.json()['result']:
            res['index'] = index_count
            index_count += 1
            res['complex_id'] = i
            print res
            #kdrew: this field makes converting to a dataframe difficult
            #res['intersections'] = ' '.join(res['intersections'])
            #kdrew: convert so that we know which genes are annotated
            res['intersection_genes'] = ' '.join([yy for j,yy in enumerate(complex_line.split()) if len(res['intersections'][j]) > 0])
            #kdrew: combine parents into a single string otherwise it gets separated into two entries, 
            #kdrew: put string back into list because pandas complains about all scalars and no index but can't set index through from_dict function
            res['parents'] = [' '.join(res['parents'])]
            del res['intersections']
            if results_df is None:
                results_df = pd.DataFrame.from_dict(res)
            else:
                results_df = pd.concat([results_df, pd.DataFrame.from_dict(res)])
    results_df = results_df.set_index("index")
    print results_df
    results_df.to_csv(args.output_filename)

        #proc = sp.Popen(['gprofiler.py', complex_line, '-c', args.correction_method, '-e', '-B', args.background_filename ], stdout=sp.PIPE, stderr=sp.PIPE)
        #gprofiler_out, err = proc.communicate()

        #print "corr. p-value\tterm ID\tt name"
        #for line in gprofiler_out.split('\n'):
        #    line_sp = line.split('\t')
        #    if len(line_sp) >= 13:
        #        if args.output_filename != None:
        #            #output_file.write("%s\t%s\t%s\n" % (line_sp[2], line_sp[8], line_sp[11]))
        #            output_file.write("%s\n" % (line,))
        #        else:
        #            #print "%s\t%s\t%s" % (line_sp[2], line_sp[8], line_sp[11])
        #            print line


    #if args.output_filename != None:
    #    output_file.close()
                


if __name__ == "__main__":
        main()

