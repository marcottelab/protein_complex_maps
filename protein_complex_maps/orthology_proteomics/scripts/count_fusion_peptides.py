import pandas as pd
import argparse



def count_peptides(inputfile, outputfile):
    fusion_full = pd.read_table(inputfile, sep="\t")
    #print fusion_full
    peps = fusion_full[['sequence']]
    #print peps

    pepgroup = peps.groupby('sequence').size()
    #for group in pepgroup:
    print pepgroup
    pepgroup.to_csv(outputfile, sep="\t", header=False)



def parse_args():

   parser=argparse.ArgumentParser(description = 'Takes an MSblender .ms2.pep.csv and counts occurences of each peptide')
   parser.add_argument('inputfile', metavar='inputfile', type=str, help = 'An MSblender output .ms2.pep.csv')
   parser.add_argument('outputfile', metavar='inputfile', type=str, help = 'An output file name')
   return parser.parse_args()


def main():
     
    args=parse_args()

    count_peptides(args.inputfile, args.outputfile)


main()
