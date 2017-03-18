import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
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

def draw_network(pairs, baitlist):
    print(pairs)
        #print(nodes)
    #https://plot.ly/ipython-notebooks/networks/
    G=nx.Graph()
   

    rowlist = []
    for index, row in pairs.iterrows(): 
         str_row = str(row['genename']) +" " + str(row['genename2']) + " " + str(row['score'])
         print(str_row)
         rowlist = rowlist + [str_row]
  
    


    G = nx.parse_edgelist(rowlist, nodetype = str, data=(('weight',float),))
    colorvalues =[]

    for node in G.nodes():
       if node in baitlist:
            colorvalues = colorvalues + ["#5ab4ac"]
            print(node, "b")
       else:
            colorvalues = colorvalues + ["#d8b365"]
            print(node, "r")
    print(colorvalues)
 
    weights = [G[u][v]['weight']*4 for u,v in G.edges()]
    pos=nx.spring_layout(G, weight=1)
    nx.draw(G, pos=pos, node_color=colorvalues, edge_color='#bdbdbd', width=weights, with_labels=True)
    plt.show()
    return



def corum_plots(data_file, inp_lower, inp_median, inp_stdev, inp_mean, lower_list, median_list, stdev_list, mean_list):

    data = np.genfromtxt(data_file,delimiter=",")
    #print(data)

    sns.set_style("whitegrid", {"axes.grid": False})
    sns.set_context("talk")
    sns.despine()

    #stdev = data[:,0]
    #median = data[:,1]
    #lower_limit = data[:,2]
    #mean = data[:,3]


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











