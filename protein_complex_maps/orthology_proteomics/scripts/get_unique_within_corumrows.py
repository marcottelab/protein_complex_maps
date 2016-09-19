#Some rows in the corum mapped to euNOGs have redundant entries
#Input line: 5919, KOG0193, KOG0193, KOG0841, KOG0841, KOG0841, KOG0841, KOG0841, KOG0841, KOG0841
#Output line: 5919,KOG0193,KOG0841

import argparse


def reducerows(infilename):


    outfilename = "nonredundant_" + infilename

    with open(infilename, "r") as infile:
        with open(outfilename, "w") as outfile:
            lines = infile.readlines()
            for line in lines:
                templine = line.strip("\n")
                linelist = templine.split(' ')
                print linelist
                uniqueline = list(set(linelist))
                finalline = " ".join(uniqueline)
                print finalline  
                outfile.write(finalline + "\n")


def main():
    parser = argparse.ArgumentParser(description='Get only unique members of each row in a space separated file')

    parser.add_argument('infilename', action="store", type=str)
    inputs = parser.parse_args()

    reducerows(inputs.infilename)




main()




