from __future__ import print_function
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import itertools as it


def main():

    parser = argparse.ArgumentParser(description="Converts a dataframe of features to libsvm format")
    parser.add_argument("--input_feature_matrix", action="store", dest="feature_matrix", required=True, 
                                    help="Filename of input feature matrix")
    parser.add_argument("--libsvm0_output_file", action="store", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")
    parser.add_argument("--libsvm1_output_file", action="store", required=False, default=None, 
                                    help="Filename of output file, default=None which prints to stdout")

    parser.add_argument("--features", action="store", dest="features", nargs='+', required=False, default=None, 
                                    help="Names of features to output, default=all")
    parser.add_argument("--label_column", action="store", dest="label_column", required=False, default='label', 
                                    help="Name of label column, default='label' if present, else 1st column")
    #parser.add_argument("--keep_labels", action="store", dest="keep_labels", nargs='+', required=False, default=None,
                                   # help="Only keep rows with these labels, default = keep all")
    parser.add_argument("--sep", action="store", dest="sep", required=False, default='$',
                                    help="Column separator for input file, default=$")

    args = parser.parse_args()

    feature_table = pd.read_csv(args.feature_matrix,sep=args.sep)

    print("opened feature matrix %s" % args.feature_matrix)
    if args.label_column in feature_table.columns:
        label_vector = feature_table[args.label_column]
        feature_table_wolabel = feature_table.drop(args.label_column,axis=1)
    else:
        #kdrew: if default label not present, default to first column 
        print("using first column as label")
        label_vector = feature_table[feature_table.columns[0]]
        feature_table_wolabel = feature_table.drop(feature_table.columns[0],axis=1)
        
    print(args.features)
    if args.features != None:
        feature_table_trim = feature_table_wolabel[args.features]
    else:
        feature_table_trim = feature_table_wolabel

    print(feature_table_trim)
    print(label_vector)


    #libsvm0_outfile = open(args.libsvm0_output_file,"wb")
    #libsvm1_outfile = open(args.libsvm1_output_file,"wb")
  

    #CDM get rows labeled one
    print("making libsvm1")
    count = 1
    one_labeled = feature_table_trim[~feature_table_trim[label]==0] 
    print("done filtering")
    print(libsvm1.head)
    for column in one_labeled.columns:
        
        if column == args.label_column:
             continue
        else:
            one_labeled["tmp"] = str(count)
            column_name = str(count) + column
            one_labeled[column_name] = one_labeled[['tmp', column]].apply(lambda x: ':'.join(x), axis=1)
            one_labeled.drop(column, 1)
            count = count + 1

    print(one_labeled)    
    #This is to move the label column to the first column
    one_labeled = one_labeled.set_index([args.label_column])
    one_labeled = one_labeled.reset_index()


    one_labeled.to_csv(args.libsvm1_output_file, index=False, header=False, sep=" ")


    #CDM get rows labeled zero
    print("making libsvm0")
    count = 1
    zero_labeled = feature_table_trim[feature_table_trim[label]==0] 
    for column in zero_labeled.columns:
        
        if column == args.label_column:
             continue
        else:
            zero_labeled["tmp"] = str(count)
            column_name = str(count) + column
            zero_labeled[column_name] = zero_labeled[['tmp', column]].apply(lambda x: ':'.join(x), axis=1)
            zero_labeled.drop(column, 1)
            count = count + 1

    
    #This is to move the label column to the first column
    print(zero_labeled)
    zero_labeled = zero_labeled.set_index([args.label_column])
    zero_labeled = zero_labeled.reset_index()
       
    zero_labeled.to_csv(args.libsvm1_output_file, index=False, header=False, sep=" ")




    #for i,index in enumerate(zero_labeled.index):
    #    libsvm0_outfile.write("%s " % label_vector[i])
    #    [libsvm0_outfile.write("%s:%s " % (i+1,x)) for i,x in enumerate(feature_table_trim.ix[index])]
    #    libsvm0_outfile.write("\n")

    #libsvm0_outfile.close()
   



    #CDM get rows labeled 1 or -1
    #one_labeled = feature_table_trim[!feature_table_trim[label]==0] 
    #for i,index in enumerate(one_labeled.index):
    #        libsvm1_outfile.write("%s " % label_vector[i])
    #        [libsvm1_outfile.write("%s:%s " % (i+1,x)) for i,x in enumerate(feature_table_trim.ix[index])]
    #        libsvm1_outfile.write("\n")

    #libsvm1_outfile.close()



    #for i,index in enumerate(feature_table_trim.index):

        #kdrew: allows the ability to only keep rows with certain labels
        #if args.keep_labels != None and str(label_vector[i]) not in args.keep_labels:
        #    continue 

     #   if label_vector[i] == 0:
      #      libsvm0_outfile.write("%s " % label_vector[i])
       #     [libsvm0_outfile.write("%s:%s " % (i+1,x)) for i,x in enumerate(feature_table_trim.ix[index])]
        #    libsvm0_outfile.write("\n")

      #  elif label_vector[i] in [1,-1]:
       #     libsvm1_outfile.write("%s " % label_vector[i])
        #    [libsvm1_outfile.write("%s:%s " % (i+1,x)) for i,x in enumerate(feature_table_trim.ix[index])]
         #   libsvm1_outfile.write("\n")


       



if __name__ == "__main__":
    main()


