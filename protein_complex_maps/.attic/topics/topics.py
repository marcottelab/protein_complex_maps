from __future__ import print_function
import argparse
import numpy as np
import lda
import lda.datasets
import pandas as pd

def format_data(infile, separator):

    df = pd.DataFrame(pd.read_table(infile, sep=separator, index_col=0))
    print(df)

    
    ids = list(df.index)
    print(ids)
    values = df.values
    print(values)

    t_values = np.transpose(values)
    print(t_values)
    print(len(ids))

    print(values.shape)
    print(t_values.shape)   
    t_values = t_values.astype(int) 
    print(t_values)
    return ids, t_values


def find_topics(ids, values, n, outfilename):
  
     print("initialize model")
     model = lda.LDA(n_topics=n, n_iter=5, random_state=1)
     print("fit model")
     model.fit(values)

     print("something")
     topic_word=model.topic_word_
     #n_top_words = 8

     with open(outfilename, "w") as outfile:
         #outfile.write(
         for i, topic_dist in enumerate(topic_word):
             #topic_words = np.array(ids)[np.argsort(topic_dist)][:-(n_top_words+1):-1]

             print(i, topic_dist)
             topic_words = np.array(ids)[np.argsort(topic_dist)]
        
             outfile.write('{}:{}{}'.format(topic_dist[0], ' '.join(topic_words), "\n"))






def main():
     parser = argparse.ArgumentParser(description="Discover protein complexes/topics")
     parser.add_argument("--infile", action="store" , required=True,
                                    help="elution profile. col1 = IDs, header = fraction")
     parser.add_argument("--sep", action="store", dest='separator', required=True,
                                     help="separator for the input elution profile")
     args = parser.parse_args()


     ids, values = format_data(args.infile, args.separator)
     for n in [10000, 15000, 20000, 25000,30000]:

         outfilename = "n" + str(n) + "_topics_" + args.infile
         find_topics(ids, values, n, outfilename)
         






main()





























