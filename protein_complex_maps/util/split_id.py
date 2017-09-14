from __future__ import print_function
import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser(description="Split an ID stored in one column to two columns")
    parser.add_argument("--input_labeled", action="store", dest="input_labeled", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default=',',
                                    help="Separator for reading feature csv, default=,")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--label_column", action="store", dest="label_column",  required=False,
                                    help="Column containing label")
    parser.add_argument("--id_column", action="store", dest="id_column",  required=False, default='ID',
                                    help="Column containing ID pairs")
    parser.add_argument("--id_sep", action="store", dest="id_sep", required=False, default=' ',
                                    help="If IDs are stored in one column, separator to split on")
    parser.add_argument("--id_columns", action="store", dest="id_columns", nargs='+', required=False, default=['key','key2'],
                                    help="List of columns that specify ids in feature matrix, default: key1 key2")

    args = parser.parse_args()

    labeled = pd.read_csv(args.input_labeled,sep=args.sep)

    only_id = labeled[args.id_column]  
    ID = args.id_column 

    #cdm Split 'IDsepID' to 'ID' 'ID'
    new_labels= only_id.str.rsplit(n=1, expand=True)
 
    new_labels.columns=['key1', 'key2']


    if args.label_column:
           only_label = labeled[args.label_column]

           new_labels = pd.concat([new_labels, only_label], axis=1)
    
    new_labels.to_csv(args.output_file,sep=" ",index=False,header=True)



if __name__ == "__main__":
    main()


