import pandas as pd
import argparse
import time
from bokeh.charts import Line, show, output_file, vplot, hplot
from bokeh.charts import defaults
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import HoverTool 
from bokeh.embed import components
import feather
import numpy as np
defaults.width = 500
defaults.height = 100


def normalize(df):
    #print "setting multiindex"
    #df = df.set_index(['fraction', 'fraction_order', 'spec', 'experiment'])
 
    result = df.copy()
    for feature_name in df.columns:
        
        max_value = df[feature_name].max()
        min_value = df[feature_name].min()
        result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
    result = result.fillna(0)
    return result


def group_normalize(df):

    df = df.set_index(['fraction', 'fraction_order', 'spec', 'experiment'])
    result = df.groupby(level = 'experiment').apply(normalize)
    print result
    return result
def tmp():

    print "setting multiindex"
    df = df.set_index(['fraction', 'fraction_order', 'spec', 'experiment'])

    result = pd.DataFrame()
     
    print "grouping"
    grouped = df.groupby(level = 'experiment')
    for name, group in grouped:
       tmp_result = pd.DataFrame()
       result_tmp = normalize(group)
       print name
       #for feature_name in group.columns:
              #group[feature_name].apply(lambda x: (x - x.mean()) / x.std())
        #       group[feature_name].apply(lambda x: (x - np.mean(x)) / (np.max(x) - np.min(x)))
         #      print "bla"
          #     print group[feature_name]
       result = result.append(result_tmp)
       #print result
    return result
    

#The other plots are linked to this first one for rescaling and zooming and such

def format_for_sparklines(elution_file):

   
    df  = pd.read_csv(elution_file)
    #df  = pd.read_csv("test_plants.csv")

   
    df = df.set_index(["ID"])
    
    df = df.transpose()

    df.to_csv("tmp")    
    frac_descs = pd.read_csv("Fraction_Details.csv")
    frac_descs  = frac_descs.set_index(['fraction'])

    frac_order = pd.read_csv("Fraction_Order.csv")
    frac_order  = frac_order.set_index(['experiment'])

    df = df.join(frac_descs)
    df.index.rename('fraction', inplace=True) 
   
    #print df
    
 
    df = df.reset_index()
    df = df.set_index(['experiment'])
    df = df.join(frac_order)
    df = df.reset_index()
       
    
    print(df)
    df = group_normalize(df)
    df = df.reset_index()


    print(df)
    
    start = time.time()
    formatted_outfile = "formatted_" + elution_file
    df.to_csv(formatted_outfile, index=False)
    end = time.time()
    print(end-start)   
    start = time.time()
    feather_outfile = formatted_outfile + ".feather"
    feather.write_dataframe(df, feather_outfile)

    end = time.time()
    print(end-start)   
 
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Feather and elution profile")
    parser.add_argument("--elution_file", action="store", dest="elution_file", required=True,
                                    help="Elution profile wide csv")

    args = parser.parse_args()

    format_for_sparklines(args.elution_file)
