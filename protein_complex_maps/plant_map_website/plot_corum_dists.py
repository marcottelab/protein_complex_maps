import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd

from bokeh.embed import components
from bokeh.plotting import figure, gridplot, GridSpec, hplot
from bokeh.charts import Histogram, output_file, show

#from bokeh.plotting import figure, 
#from bokeh.resources import INLINE
#from bokeh.util.string import encode_utf8
from bokeh.models import HoverTool, ColumnDataSource, Span
#from bokeh.util.browser import view




def interaction_share(data_dict1, title1, data_dict2, title2):

    sns.set_style("whitegrid", {"axes.grid": False})
    sns.set_context("talk")

    plt.subplot(2,1,1)
    plt.bar(range(len(data_dict1)), data_dict1.values(), align='center')
    plt.xticks(range(len(data_dict1)), data_dict1.keys())
    plt.title(title1)
    
    plt.subplot(2,1,2)
    plt.bar(range(len(data_dict2)), data_dict2.values(), align='center')
    plt.xticks(range(len(data_dict2)), data_dict2.keys())
    plt.title(title2)
 

    plt.show()

def get_network(pairs, id1, id2):
    #print(pairs)
        #print(nodes)
    
    #https://plot.ly/ipython-notebooks/networks/
    G=nx.Graph()
   

    rowlist = []



    for index, row in pairs.iterrows(): 
         if str(row[id1]) == "":
               val1 = "None"
         else:
               val1 = str(row[id1])
         if str(row[id2]) == "":
               val2 = "None"
         else:
               val2 = str(row[id2])


         str_row = val1 +" " + val2 + " " + str(row['score'])
         #print(str_row)
         rowlist = rowlist + [str_row]
    #print(rowlist)  
    


    G = nx.parse_edgelist(rowlist, nodetype = str, data=(('weight',float),))
    #print(dir(G))

    return G


def filter_nodes(G, mindegree = 1):
    #print(colorvalues)
    #print(G.edges)
    print "degree troubleshooting"
    print mindegree
    print type(mindegree)
    outdeg = G.degree()
    print G.degree()
    to_remove = [n for n in outdeg if outdeg[n] < int(mindegree)]
    G.remove_nodes_from(to_remove)
    print "after removing"
    print G.degree()
    print G.nodes()

    return G

def draw_network(G, baitlist, showplot=False):
   
    colorvalues =[]

    for node in G.nodes():
       if node in baitlist:
            colorvalues = colorvalues + ["#5ab4ac"]
            #print(node, "b")
       else:
            colorvalues = colorvalues + ["#d8b365"]
            #print(node, "r")




    weights = [G[u][v]['weight']*4 for u,v in G.edges()]
    pos=nx.spring_layout(G, weight=1)
    nx.draw(G, pos=pos, node_color=colorvalues, edge_color='#bdbdbd', width=weights, with_labels=True)
    if showplot == True:
        plt.show()
    return weights, pos, colorvalues





def bokeh_network(G, weights, pos, colorvalues, annot_dict):
        tools =   'wheel_zoom,reset, hover, pan, save'
        labels = [ str(v) for v in G.nodes() ]
        annotations = []
        for lab in labels:
             annotations.append(annot_dict[lab])
        print annotations


        print labels
        vx, vy = zip(*[ pos[v] for v in G.nodes() ])
        xs, ys = [], []
        for (a, b) in G.edges():
            x0, y0 = pos[a]
            x1, y1 = pos[b]
            xs.append([x0, x1])
            ys.append([y0, y1])
        print ("positions")
        print pos
        print xs, ys 
        hover = HoverTool(names=["circles"])

        f = figure(plot_width=800, plot_height=700,
                   x_axis_type=None, y_axis_type=None,
                   outline_line_color=None,
                   tools=tools, x_range=[-0.1, 1.1])

        sourceDS = ColumnDataSource(
            data=dict(
                xname=vx,
                yname=vy,
                name = labels,
                annot = annotations
             )
        )
        print sourceDS
        hover= f.select_one(HoverTool)
        #hover.tooltips = [
        #("name", "@name"), ("degree", "@degree_label"), ("annotation", "@annot")]
        hover.tooltips = [
        ("Annotation", "@annot")]


         #f.circle(vx, vy, size=18, name="circles",line_color=colorvalues, fill_color=colorvalues)
        f.circle('xname', 'yname', source=sourceDS, name="circles", size=21 ,line_color=colorvalues, fill_color=colorvalues)
        f.multi_line(xs, ys, line_color="#bdbdbd", line_width=weights)
        f.text(vx, vy, text=labels,
               text_font_size="18px", text_align="center", text_baseline="middle")


        return components(f)

def hist_and_spans(data, stat, input_stat, stat_list, tools):
   
    data = pd.DataFrame(data[stat].dropna())
    #data = data[stat].replace('NA', 'NaN')
    #data = data.fillna('')
    print data.to_string()
    hist = figure(width=300, height=300, tools = tools)
    hist = Histogram(data, values=stat, bins=50)
    #Add vlines inp_median
    for s in stat_list: 
        if s:
            hist.add_layout(Span(location=s, line_color='green', dimension='height'))
    if input_stat:    
        hist.add_layout(Span(location=input_stat, line_color='red', dimension='height'))


    return hist
def corum_plots_bokeh(data_file,  inp_lower, inp_median, inp_stdev, inp_mean, lower_list, median_list, stdev_list, mean_list):
    tools = 'xwheel_zoom,reset, hover'
    #data = np.genfromtxt(data_file,delimiter=" ")
    data = pd.DataFrame(pd.read_csv(data_file, sep=",", header=None)) 
    data.columns = ['lower', 'median', 'stdev', 'mean']
    print(data)

    hist1 = hist_and_spans(data, 'median', inp_median, median_list, tools)
    hist2 = hist_and_spans(data, 'mean', inp_mean, mean_list, tools)
    print inp_stdev, stdev_list
    print inp_lower, lower_list
   
    #hist3 = hist_and_spans(data, 'stdev', inp_stdev, stdev_list, tools)
    #hist4 = hist_and_spans(data, 'lower', inp_lower, lower_list, tools)
    p = hplot(hist1, hist2)
    #p = gridplot([hist1, hist2], [hist3, hist4])
   #plt.axvline(x=inp_median, color='r')
    #for med in median_list:
    #   plt.axvline(x=med, color='r', alpha=0.4)
    return components(p)
 #   output_file("histogram_single.html", title="histogram_single.py example")

#    show(hist4)



 
def corum_plots(data_file, inp_lower, inp_median, inp_stdev, inp_mean, lower_list, median_list, stdev_list, mean_list):
    #print(data_file)
    data = np.genfromtxt(data_file,delimiter=" ")
    #print(data)

    sns.set_style("whitegrid", {"axes.grid": False})
    sns.set_context("talk")
    sns.despine()

    #stdev = data[:,0]
    #median = data[:,1]
    #lower_limit = data[:,2]
    #mean = data[:,3]
    print("data for plotting")
    #print(data)

    num_bins = 20
    # the histogram of the data

    plt.subplot(2, 2, 1)
    plt.hist(data[:,1], num_bins)
    plt.xlabel('Median_Score')
    plt.ylabel('Num Complexes')
    plt.xlim(0, 1)
    plt.axvline(x=inp_median, color='r')
    for med in median_list:
       plt.axvline(x=med, color='r', alpha=0.4)
  

    sns.despine()  

    plt.subplot(2, 2, 2)   
    plt.hist(data[:,3], num_bins)
    plt.xlabel('Mean_Score')
    plt.ylabel('Num Complexes')
    plt.xlim(0, 1)
    plt.axvline(x=inp_mean, color='r')
    for m in mean_list:
       plt.axvline(x=m, color='r', alpha=0.4)
 

    sns.despine()  

    plt.subplot(2, 2, 3)   
    plt.hist(data[:,0], num_bins)
    plt.xlabel('Lower_Limit: Median - StDev')
    plt.ylabel('Num Complexes')
    plt.xlim(0, 1)
    plt.axvline(x=inp_lower, color='r')
    for low in lower_list:
       plt.axvline(x=low, color='r', alpha=0.4)
 
    sns.despine()  

    plt.subplot(2, 2, 4)   
    plt.hist(data[:,2], num_bins)
    plt.xlabel('Standard_Deviation')
    plt.ylabel('Num Complexes')
    plt.xlim(0, 1)
    plt.axvline(x=inp_stdev, color='r')
    for sd in stdev_list:
       plt.axvline(x=sd, color='r', alpha=0.4)
 
    sns.despine()  
    
    # Tweak spacing to prevent clipping of ylabel
    #plt.subplots_adjust(left=0.15)
    plt.show()
    


def load_args():
 

    parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
    parser.add_argument("--data", action="store", dest="data", required=True,
                  help="Four comma separated columns stdev,median,lower_limit,mean")

    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = load_args()
    
    corum_plots(args.data)











