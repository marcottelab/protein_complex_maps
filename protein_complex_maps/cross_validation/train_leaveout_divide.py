from __future__ import print_function
import random
import argparse

<<<<<<< HEAD

def partition ( lst, n ):
    #From def partition ( lst, n ):
    return [ lst[i::n] for i in xrange(n) ]


def make_divide(row_file, n, seed):

    random.seed(seed)

    with open(row_file, "r") as handle:
        row_list = handle.readlines()
        random.shuffle(row_list) #This is inplace
        print(len(row_list))
        divisions = partition(row_list, n)
        print(len(divisions))
        for i in range(0,n):
            test = divisions[i]
            train= [] 
            for j in divisions:
                if j != test:
                    train = train + j

            outfilename_test = "_".join(["leaveout_seg", str(i+1), row_file])

            outfile_test = open(outfilename_test, "w")
            for label in test:
                 outfile_test.write(label)

            outfilename_train = "_".join(["train_seg", str(i+1), row_file])

            outfile_train = open(outfilename_train, "w")
            for label in train:
                 outfile_train.write(label)

            outfile_test.close()
            outfile_train.close()


            #train = [j for j in divisions if j!=i]
            #train = [item for sublist in train for item in sublist]
            print(len(test))
            print(len(train))             


def load_file():

=======
import pandas as pd

def divide_file(input_df, n, seed, label="label"):
    #kdrew: get only positive and negative labeled rows
    df_pos = input_df.query("%s == 1" % label)
    df_neg = input_df.query("%s == -1" % label)

    #kdrew: shuffle positive and negative dataframe and get indices
    pos_rand_ids = df_pos.sample(frac=1, random_state=seed).index
    neg_rand_ids = df_neg.sample(frac=1, random_state=seed).index

    #kdrew: split indices into separate partitions
    pos_rand_id_splits = [ pos_rand_ids[i::n] for i in xrange(n) ]
    neg_rand_id_splits = [ neg_rand_ids[i::n] for i in xrange(n) ]

    #kdrew: combine indices of all other sets, leaving out current set
    pos_combined_splits = []
    neg_combined_splits = []
    for i in xrange(n):
        pos_combined = []
        neg_combined = []
        for j in xrange(n):
            if i != j:
                pos_combined = pos_rand_id_splits[j].union(pos_combined)
                neg_combined = neg_rand_id_splits[j].union(neg_combined)
        pos_combined_splits.append(pos_combined)
        neg_combined_splits.append(neg_combined)
    
    #kdrew: store split dataframes
    df_pos_training = [df_pos.ix[pos_combined_splits[i]] for i in xrange(n)]
    df_neg_training = [df_neg.ix[neg_combined_splits[i]] for i in xrange(n)]
    df_pos_leaveout = [df_pos.ix[pos_rand_id_splits[i]] for i in xrange(n)]
    df_neg_leaveout = [df_neg.ix[neg_rand_id_splits[i]] for i in xrange(n)]

    return_dict = dict()
    return_dict['pos_train'] = df_pos_training
    return_dict['neg_train'] = df_neg_training
    return_dict['pos_leaveout'] = df_pos_leaveout
    return_dict['neg_leaveout'] = df_neg_leaveout
    return return_dict


def main():
>>>>>>> 28fde78333cd780321f18246dc56cd4ef64322d7

    parser = argparse.ArgumentParser(description='Divide files into test and train set')
    parser.add_argument('--input_file', action="store", type=str, required=True,
                            help="File of rows to be divided")
    parser.add_argument('--output_filename', action="store", type=str, required=True,
                            help="Output filename, string gets mangled to generate multiple output files")
    parser.add_argument('--sep', action="store", type=str, default=',', required=False,
                            help="Separator of input file")
    parser.add_argument('--n', action="store", type=int,  default=5, required=False, help="Number of chunks to divide file into")
    parser.add_argument('--seed', action="store", type=int,  default=42, required=False, help="Random seed")
    parser.add_argument("--id_columns", action="store", nargs='+', dest="id_columns", required=True, 
                                    help="columns that specifies the ids in feature matrix")

    args = parser.parse_args()
    #make_divide(args.input_file, args.n, args.seed)
    input_df = pd.read_csv(args.input_file, sep=args.sep)
    divided_file_dict = divide_file(input_df, args.n, args.seed)

    #kdrew: output each training and leaveout set for both postives and negatives
    for i in xrange(args.n):
        #kdrew: output concatenated positive and negative training
        df_pos = divided_file_dict['pos_train'][i]
        df_neg = divided_file_dict['neg_train'][i]
        df = pd.concat([df_pos, df_neg])
        specific_output_filename = args.output_filename + ".train" + str(i) + ".txt"
        df.to_csv(specific_output_filename, sep=args.sep, index=False)

        #kdrew: output concatenated positive and negative leaveout
        df_pos = divided_file_dict['pos_leaveout'][i]
        df_neg = divided_file_dict['neg_leaveout'][i]
        df = pd.concat([df_pos, df_neg])
        specific_output_filename = args.output_filename + ".leaveout" + str(i) + ".txt"
        df.to_csv(specific_output_filename, sep=args.sep, index=False)

    #kdrew: output leaveout ppis
    for i, df in enumerate(divided_file_dict['pos_leaveout']):
        specific_output_filename = args.output_filename + ".pos_leaveout_ppis%s.txt" % i
        df[args.id_columns].to_csv(specific_output_filename, sep=' ', index=False, header=False)
    for i, df in enumerate(divided_file_dict['neg_leaveout']):
        specific_output_filename = args.output_filename + ".neg_leaveout_ppis%s.txt" % i
        df[args.id_columns].to_csv(specific_output_filename, sep=' ', index=False, header=False)



if __name__ == "__main__":
    main()
