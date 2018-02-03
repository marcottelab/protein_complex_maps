from __future__ import print_function
import argparse
import pandas as pd


def annotate_IDs(pairs, annots):

    pairs=pairs.rename(columns = {0:'ID1'})
    pairs=pairs.rename(columns = {1:'ID2'})
  
    pairs = pd.merge(pairs, annots, left_on='ID1', right_on='ID', how='left')    
    pairs = pd.merge(pairs, annots, left_on='ID2', right_on='ID', how='left', suffixes=['_1','_2'])    
   
    pairs=pairs.drop(['ID_1', 'ID_2'], axis=1)
 
    return pairs

def main():

    parser = argparse.ArgumentParser(description="With IDs in first two columns, retrieve annotations for each ID")
    parser.add_argument("--input_pairs", action="store", dest="input_pairs", required=True, 
                                    help="Filename of input pairs, (first two columns IDs) No header")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=' ',
                                    help="Separator for reading feature csv, default=$")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--input_annots", action="store", dest="input_annots", required=True, 
                                    help="Filename of annotation file")
    parser.add_argument("--annot_sep", action="store", dest="annot_sep", required=False, default=',', 
                                    help="Sep of annotation file")
    parser.add_argument("--annot_key", action="store", dest="annot_key", required=False, default='ID', 
                                    help="Col name in annotation file that matches first two rows of input file")
    parser.add_argument('-t','--target_col', action="store", dest="target_col",  required=False, default=None, help="designate a specific column to full from annotation file")

    parser.add_argument("--output_sep", action="store", dest="output_sep", required=False, default=' ', 
                                    help="Sep of output file")
  
    args = parser.parse_args()

    pairs = pd.read_csv(args.input_pairs,sep=args.sep, header=None)
    annots = pd.read_csv(args.input_annots,sep=args.annot_sep)

    if args.target_col:
       annots = annots[[args.annot_key, args.target_col]]

    annotated = annotate_IDs(pairs,annots)
 
    annotated.to_csv(args.output_file,sep='\t',index=False,header=True)

if __name__ == "__main__":
    main()


