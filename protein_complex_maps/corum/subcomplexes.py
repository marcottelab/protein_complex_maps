from __future__ import print_function
import csv
import argparse    

def get_subcomplexes(infilename, outfilename, sep):

    complexes_list = []
    complex_dict = {}
    

    with open(infilename,'rb') as infile:
        complexes = csv.reader(infile, delimiter=sep)
        n = 1
        for row in complexes:
            complexes_list.append(row)
            complex_dict[n] = row 
            n = n + 1
    #print complex_dict
    #print(complexes_list)  
    subcomplex_dict = {}
    to_remove = []
    for complex1 in complexes_list:
        for complex2 in complexes_list:
            if complex1 == complex2:
               continue
            elif all(x in complex1 for x in complex2):
                print("%s : %s\n" % (complex2, complex1))
                to_remove.append(complex2)
    to_remove = [list(x) for x in set(tuple(x) for x in to_remove)]
    
    print(to_remove)
    final_list = []
    with open(outfilename, "w") as outfile:
        for c in complexes_list:
            if c not in to_remove:
                final_list.append(c)
                final_string = ' '.join(c) + '\n'
                print(final_string)
                outfile.write(final_string)
 
    print("%s complexes in original list" % str(len(complexes_list)))
    print("%s complexes removed" % str(len(to_remove)))
    print("%s final complexes" % str(len(final_list)))
       
 

parser = argparse.ArgumentParser(description="")

parser.add_argument('--corum_file', action="store", type=str, help="one group per line")
parser.add_argument('--outfile', action="store", type=str, help="outfile name")
parser.add_argument('--sep', action="store", type=str, default=' ', required=False, help="separator for input and output file")
inputs = parser.parse_args()

get_subcomplexes(inputs.corum_file, inputs.outfile, inputs.sep)



