#! /usr/bin/env python
from __future__ import print_function
import argparse
from protein_complex_maps.features.ExtractFeatures.Features import Elut

parser = argparse.ArgumentParser(description="Print information on one or several elution experiment(s)")
parser.add_argument("infiles", nargs="+",help='Elution profile csv file. Optionally takes a wildcard like *.csv')
args = parser.parse_args()

def pretty_dict(D,indent_width=0):
    for key,val in D.iteritems():
        if type(val) is dict:
            print(' '*indent_width + key)
            pretty_dict(val,indent_width*2)
        else:
            print(' '*indent_width + "\t".join([key,str(val)]))
 
if __name__ == '__main__':
    for f in args.infiles:
        e = Elut()
        e.load(f)
        print(f)
        pretty_dict(e.info,4)
        if len(args.infiles) > 1:
            print
    
