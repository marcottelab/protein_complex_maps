import os.path
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

def main():

    parser = argparse.ArgumentParser(description="Tool to read in msblender file and produce elution file")
    parser.add_argument("--prot_count_files", action="store", dest="prot_count_files", nargs='+', required=True, 
                        help="Filenames of MSblender prot_count files")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True,
                        help="Filenames of MSblender prot_count files")
    parser.add_argument("--spectral_count_type", action="store", dest="spectral_count_type", required=False, default='Unique',
                        help="Which spectral count type should be stored? Unweighted, Weighted or Unique (default)")
    parser.add_argument("--remove_zero_unique", action="store_true", dest="remove_zero_unique", required=False, default=False,
                        help="Flag for removing entries with zero unique matches, default=False")
    parser.add_argument("--fraction_name_from_filename", action="store_true", dest="fraction_name_from_filename", required=False, default=False,
                        help="Flag for parsing fraction name from filename, default=False, defaults to sequential numbering")
    
    # Formatting flags
    parser.add_argument("--contaminant_flag", action='store', dest='contaminant_flag', required=False, default="CONTAMINANT",
                        help="If specified, remove rows with this string from output file. Can specify None if you want to retain all rows for some reason, default=CONTAMINANT")
    parser.add_argument("--msblender_format", action="store_true", dest="msblender_format", required=False, default=False,
                        help="DEPRECATED: Original file format has the first column as the sum of all the columns, default=False")
    parser.add_argument("--parse_uniprot_id", action="store_true", dest="parse_uniprot_id", required=False, default=False,
                        help="Parse out the uniprot id from the index and store in index, default=False")

    args = parser.parse_args()

    elution_profile_df = pd.DataFrame() 
    
    for i, filename in enumerate(args.prot_count_files):
        print "Reading {}".format(filename)
        df = pd.read_table(filename, index_col=0)
        
        ## BJL: Useful code if we want to eventually parse other MSBlender outfiles besided .group
        #else:
        #   #kdrew: the non-group files have an extra line at the bottom for total counts
        #   df = pd.read_table(filename, index_col=0, skip_footer=1)
        
        single_elution_df = pd.DataFrame(df[args.spectral_count_type])

        if args.remove_zero_unique:
            single_elution_df = single_elution_df[single_elution_df['Unique'] > 0]

        if args.fraction_name_from_filename:
            fraction_name = os.path.basename(filename).split('.')[0]
            #print fraction_name # Debugging

            single_elution_df.columns = [fraction_name]
        else:
            single_elution_df.columns = [i]

        #print single_elution_df # Debugging
        elution_profile_df = pd.concat([elution_profile_df, single_elution_df], axis=1)


    elution_profile_df = elution_profile_df.fillna(0)
    
    if args.contaminant_flag != None:
        elution_profile_df = elution_profile_df[~elution_profile_df.index.str.contains(args.contaminant_flag)]

    if args.parse_uniprot_id:
        uniprots = [x.split('|')[1] if len(x.split('|'))>1 else x for x in elution_profile_df.index ]
        elution_profile_df.index = uniprots

    if args.msblender_format:
        print "WARNING: --msblender_format has been deprecated to put Blake's scripts on an iron lung"
        total_count = pd.DataFrame(elution_profile_df.sum(axis=1))
        total_count.columns = ['TotalCount']
        elution_profile_df = pd.concat([total_count, elution_profile_df], axis=1)
        elution_profile_df.index.name='#ProtID'
        print elution_profile_df.index.name

    #print elution_profile_df.index.name # Debugging
    elution_profile_df.to_csv(args.output_filename, sep='\t')


if __name__ == "__main__":
    main()



