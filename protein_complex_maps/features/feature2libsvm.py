
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it


def main():

    parser = argparse.ArgumentParser(description="Converts a dataframe of features to libsvm format")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--output_file", action="store", dest="out_filename", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--features", action="store", dest="features", nargs='+', required=False, default=None, 
                                    help="Names of features to output, default=all")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else 1st column")
    parser.add_argument("--keep_labels", action="store", dest="keep_labels", nargs='+', required=False, default=None,
                                    help="Only keep rows with these labels, default = keep all")

    args = parser.parse_args()


    feature_table = pd.read_csv(args.feature_matrix,sep='$')

    if args.label_column in feature_table.columns:
        label_vector = feature_table[args.label_column]
        feature_table_wolabel = feature_table.drop(args.label_column,axis=1)
    else:
        #kdrew: if default label not present, default to first column 
        print "using first column as label"
        label_vector = feature_table[feature_table.columns[0]]
        feature_table_wolabel = feature_table.drop(feature_table.columns[0],axis=1)
        

    if args.features != None:
        feature_table_trim = feature_table_wolabel[args.features]
    else:
        feature_table_trim = feature_table_wolabel

    print feature_table_trim
    print label_vector

    outfile = open(args.out_filename,"wb")

    for i,index in enumerate(feature_table_trim.index):

        #kdrew: allows the ability to only keep rows with certain labels
        if args.keep_labels != None and str(label_vector[i]) not in args.keep_labels:
            continue 

        outfile.write("%s " % label_vector[i])
        #print feature_table_trim.ix[index]
        [outfile.write("%s:%s " % (i+1,x)) for i,x in enumerate(feature_table_trim.ix[index])]
        outfile.write("\n")

    outfile.close()
        



if __name__ == "__main__":
    main()


