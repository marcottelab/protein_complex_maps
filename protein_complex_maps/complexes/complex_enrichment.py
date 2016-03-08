    
import numpy as np
import argparse
import itertools as it

import subprocess as sp



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

    args = parser.parse_args()

    #ids = "Q13772 Q13200 O00487 Q9BRP4 O43242 P62195 P43686 Q99460 P62333 P35998 P62191 Q16186 P17980 Q05086 Q15008 Q9Y5K5 Q9UNM6 O00231 O00232 Q53HC0 O75832 P51665 P48556 P55036 Q16401 O00233 Q96EN9"

    if args.output_filename != None:
        output_file = open(args.output_filename,'wb')

    complex_file = open(args.complex_filename,'rb')
    for i, complex_line in enumerate(complex_file.readlines()):
        if args.output_filename != None:
            output_file.write("complex: %s\n" % (i,))
            output_file.write("""#  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list\n""")
        else:
            print "complex: %s" % (i,)
            print """#  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list"""
        
        proc = sp.Popen(['gprofiler.py', complex_line, '-c', 'bonferroni', '-e', '-B', args.background_filename ], stdout=sp.PIPE, stderr=sp.PIPE)
        gprofiler_out, err = proc.communicate()

        #print "corr. p-value\tterm ID\tt name"
        for line in gprofiler_out.split('\n'):
            line_sp = line.split('\t')
            if len(line_sp) >= 13:
                if args.output_filename != None:
                    #output_file.write("%s\t%s\t%s\n" % (line_sp[2], line_sp[8], line_sp[11]))
                    output_file.write("%s\n" % (line,))
                else:
                    #print "%s\t%s\t%s" % (line_sp[2], line_sp[8], line_sp[11])
                    print line


    output_file.close()
                


if __name__ == "__main__":
        main()

