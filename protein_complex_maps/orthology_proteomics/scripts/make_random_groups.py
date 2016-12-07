from __future__ import print_function
import argparse
import random

def make_random_corum(infile, outfile, id_pool, sep):
   ''' This function makes random complexes that match the size
       of a given list of corum complexes.

   '''

   #output from get_elution_ids.py

   out = open(outfile, "w")
   pool = open(id_pool, "r").readlines()


   for raw_line in open(infile, "r").readlines():
       line = raw_line.split(sep)
       print(line)
       complex_size = len(line)
       
       rand_complex = random.sample(pool, complex_size)
       print(rand_complex)
       outline = ' '.join(rand_complex)
       outline = outline.replace("\n", "")
       out.write(outline + "\n")       
 
     
   out.close()
   #pool.close()
       
 
parser = argparse.ArgumentParser(description='Make a bunch of random complexes to match corum')

parser.add_argument('--corum_file', action="store", type=str, help="one group per line")
parser.add_argument('--outfile', action="store", type=str, help="random groups drawn from id_pool")
parser.add_argument('--id_pool', action="store", type=str, help="One Identifier per line")
parser.add_argument('--sep', action="store", type=str, default=' ', required=False, help="separated for reading and writing files")
inputs = parser.parse_args()

make_random_corum(inputs.corum_file, inputs.outfile, inputs.id_pool, inputs.sep)






    
