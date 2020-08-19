
import numpy as np
import argparse
import pandas as pd
import pandas as pd #for averaging dictionary values
import matplotlib.pyplot as plt 

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import average_precision_score



def main():

    parser = argparse.ArgumentParser(description="Generate Precision Recall Curve")
    parser.add_argument("--results_wprob", action="store", dest="results_wprob", nargs='+', required=True, 
                                    help="Filename of pair results with probability")
    parser.add_argument("--input_positives", action="store", dest="positives", nargs='+', required=True, 
                                    help="Filename of positive pairs")
    parser.add_argument("--input_negatives", action="store", dest="negatives", nargs='+',
                                    help="Filename of negative pairs, default = None (generated from processing positives)")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--labels", action="store", dest="labels", nargs='+', required=False, default=None,
                                    help="Labels for input results in order of sets of pairs")
    parser.add_argument("--threshold", action="store", type=float, dest="threshold", required=False, 
                                    help="Only tally predictions above probability threshold")
    parser.add_argument("--plot_thresholds", action="store_true", dest="plot_thresholds", required=False, default=False,
                                    help="Add probability threshold markers to plot")
    parser.add_argument("--average_probs", action="store_true", dest="avg_probs", required=False, default=False,
                                    help="Average a set of result_wprob instead of plotting individually")

    parser.add_argument("--header", action="store_true", dest="header", required=False, default=False,
                                    help="Set true in results_wprob has a header")
    parser.add_argument("--id_cols", action="store", nargs = "+", dest="id_cols", required=False, type = int, default=[0,1],
                                    help="Which column(s) in input results are the identifiers")
    parser.add_argument("--prob_col", action="store", dest="prob_col", required=False, type =int,default=2,
                                    help="Which column in input results contains the probability score")
    parser.add_argument("--results_delim", action = "store", required = False, default = "\t",                                                                                    help="separator for input results")
    parser.add_argument("--ppi_sep", action="store", required=False,
                                    help="If ppi identifiers are separated")
    parser.add_argument("--labels_delim", action = "store", required = False, default = "\t",                                                                                    help="separator for input positives and negatives")
    parser.add_argument("--plot",  dest="plot", required=False, default=False,choices=('True','False'),
                                    help="Plot pr plot")
 

    args = parser.parse_args()

    # Not using bool due to argparse handling bools
    if args.plot == 'True':
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.mlab as mlab
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        rcParams.update({'figure.autolayout': True})
        mpl.rc('pdf', fonttype=42)
        import seaborn as sns
        sns.set_style("white")
        

    output_pr_table = pd.DataFrame()
    #output_roc_table = pd.DataFrame()


    # Parse and organize results  
    for result_file in args.results_wprob:

        print(result_file)
        if args.header:
            results = pd.read_csv(result_file, sep = args.results_delim, engine='python') 
        else:
            results = pd.read_csv(result_file, header = None, sep = args.results_delim, engine='python')
    
        if len(args.id_cols) == 1:
          if not args.ppi_sep:
             print("--ppi_sep required for PPI stored in one column")
             return 
          else:
           results['ID'] =results.iloc[:, args.id_cols[0]]
           results['ID'] = results.apply(lambda row : frozenset(row['ID'].split(args.ppi_sep)), axis = 1)
    
    
        elif len(args.id_cols) == 2:
           results.ID1 = results.iloc[:, args.id_cols[0]]
           results.ID2 = results.iloc[:, args.id_cols[1]]
           results['ID'] = results.apply(lambda row : frozenset([row['ID1'], row['ID2']]), axis = 1)
            
    
     
       
        results['prob'] = results.iloc[:,args.prob_col]
        results = results.set_index(["ID"])
        print(results)
        for i in range(len(args.positives)):
    
           positives = pd.read_csv(args.positives[i], header = None, names = ["ID1", "ID2"], sep = args.labels_delim, engine='python')
           print(positives)
    
           positives['numlabel'] = 1
    
           negatives = pd.read_csv(args.negatives[i], header = None, names = ["ID1", "ID2"], sep = args.labels_delim, engine='python')
           negatives['numlabel'] = -1
           print(negatives)       

           labels = pd.concat([positives,negatives])
           labels['ID'] = labels.apply(lambda row : frozenset([row['ID1'], row['ID2']]), axis = 1)
    
     
           if args.labels:
              labels.label = args.labels[i]
      
           labels = labels.set_index(["ID"])
     
           labelled_results= results.join(labels,  how = "inner", rsuffix = "_2")
           precision, recall, thresholds = precision_recall_curve(labelled_results.numlabel, labelled_results.prob)
           #FPR, TPR, roc_thresholds = roc_curve(labelled_results.numlabel, labelled_results.prob)
    
   
           pr_table = pd.DataFrame({'Recall':recall, 'Precision':precision, 'threshold':np.append(thresholds, np.array([1])), 'label':args.labels[i], 'file':result_file })
           pr_table['FDR'] = 1 - pr_table['Precision']
 

           output_pr_table = output_pr_table.append(pr_table)
 
           if args.plot == 'True':
               line, = plt.plot(recall, precision, label=args.labels[i])
    if args.plot == 'True':
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall')

        plt.savefig(args.output_file)

 
    output_pr_table.to_csv(args.output_file, index = False)

  

if __name__ == "__main__":
    main()


