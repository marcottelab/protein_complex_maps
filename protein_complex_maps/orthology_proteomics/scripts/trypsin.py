'''
    this is an in silico trypsin digestion program. The input is a fasta file which contains protein sequence to be digested, the output is a txt file which contains all trypsin digested peptides and corresponding protein accessions.

From https://github.com/yafeng/trypsin/blob/master/trypsin.py
Modified by Claire McWhite 5/13/2016 to make csvs
'''

def TRYPSIN(proseq,miss_cleavage):
    peptides=[]
    cut_sites=[0]
    for i in range(0,len(proseq)-1):
        if proseq[i]=='K' and proseq[i+1]!='P':
            cut_sites.append(i+1)
        elif proseq[i]=='R' and proseq[i+1]!='P':
            cut_sites.append(i+1)
    
    if cut_sites[-1]!=len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites)>2:
        if  miss_cleavage==0:
            for j in range(0,len(cut_sites)-1):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])



        elif miss_cleavage==1:
            for j in range(0,len(cut_sites)-2):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
            
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])

        elif miss_cleavage==2:
            for j in range(0,len(cut_sites)-3):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+3]])
            
            peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
            peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    else: #there is no trypsin site in the protein sequence
        peptides.append(proseq)
    final_peptides = [i for i in peptides if len(i)>=6 and len(i)<=60]
    return final_peptides

import sys
import os
import getopt
from Bio import SeqIO

if __name__ == "__main__":
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! wrong command, please read the mannual in Readme.txt."
        print "Example: python trypsin.py --input input_filename --output output_filename --miss 1"
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                             'miss=',
                                                             'output='])
        for opt, arg in options:
            if opt == '--input': input_file=arg
            elif opt == '--miss': n=int(arg)  #number of miss cleavage allowed
            elif opt == '--output':output_file=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    
    handle=SeqIO.parse(input_file,'fasta')
    output=open(output_file,'w')
    
    output.write("%s,%s\n" % ("ProteinID","Peptide"))
    
    for record in handle:
        proseq=str(record.seq)
        peptide_list=TRYPSIN(proseq,n)
        for peptide in peptide_list:
            output.write("%s,%s\n" % (record.id,peptide))
    
    handle.close()
    output.close()
    
