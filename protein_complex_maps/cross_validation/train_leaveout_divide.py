from __future__ import print_function
import random
import argparse


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

            outfilename_test = "_".join(["test_seg", str(i+1), row_file])

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


    parser = argparse.ArgumentParser(description='Divide files into test and train set')
    parser.add_argument('--input_file', action="store", type=str, help="File of rows to be divided")
    parser.add_argument('--n', action="store", type=int,  default=5, required=False, help="Number of chunks to divide file into")
    parser.add_argument('--seed', action="store", type=int,  default=42, required=False, help="Random seed")

    inputs = parser.parse_args()
    make_divide(inputs.input_file, inputs.n, inputs.seed)


load_file()













