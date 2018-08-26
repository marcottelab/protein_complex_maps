import argparse
import pandas as pd

def main():
    
    parser = argparse.ArgumentParser(description="Tidy a single elution files. From Rows = ID & columns = fractionID to long form where each row = Fraction,ID,value")
    parser.add_argument("--input_elution", action='store', dest='input_elution', required=True, help="Input csv elution files")
    parser.add_argument("--outfile", action="store", dest="outfile", required=True, help="Name of outfile")
    parser.add_argument("--firstcol_name", action = "store", dest = "firstcol_name", required = False, default = "ID", 
                        help = "Name for output index column, default = 'ID'")
    parser.add_argument("--valuecol_name", action = "store", dest = "valuecol_name", required = False, default = "value", 
                        help = "Name for output value column, default = 'value'")
    parser.add_argument("--experiment_id", action="store", dest="experiment_id", required=False, help="If present, adds ExperimentID column")
    parser.add_argument("--keepzeros", action="store_true", dest="keepzeros", default = False, required=False, help="Keep rows for fractions with zero count of the peptide")
    parser.add_argument("--elution_sep", action='store', dest='sep', required=False, default="\t", help="Separator used for input elution profiles, default='\t'")
    
    args = parser.parse_args()
    df = pd.read_csv(args.input_elution, sep = args.sep)   
   
    #Similar to tidyr::gather
    outdf = pd.melt(df, id_vars=df.columns[0], value_vars=df.columns[1:].tolist()) 
    outdf.columns = [args.firstcol_name, 'FractionID', args.valuecol_name]
    outdf = outdf[outdf[args.valuecol_name] != 0]

    if args.experiment_id:
        outdf['ExperimentID'] = args.experiment_id
        outdf = outdf[['ExperimentID', 'FractionID', args.firstcol_name, args.valuecol_name]]
    else:
        outdf = outdf[['FractionID', args.firstcol_name, args.valuecol_name]] 

    outdf.to_csv(args.outfile, index=False)
    
if __name__ == '__main__':
    main()
