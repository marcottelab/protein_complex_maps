from __future__ import print_function
from sklearn.datasets import load_svmlight_file
from sklearn.datasets import load_boston

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
#from sklearn.ensemble import BaggingClassifier
#from sklearn.ensemble import BaggingRegressor
#from sklearn.ensemble import VotingClassifier

from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
import numpy as np
import argparse
import pandas as pd
#Load boston housing dataset as an example



def extratrees(X_dense, Y):    
    print("start extra trees")
    etc = ExtraTreesClassifier()
    etc.fit(X_dense, Y)
    # display the relative importance of each attribute
    print(etc.feature_importances_)
    print(sorted(zip(map(lambda x: round(x, 4), etc.feature_importances_), names), 
             reverse=True))



def randomforest(X_dense, Y):

    print("fitting forest")

    rf = RandomForestRegressor()
    rf.fit(X_dense, Y)
    print("Features sorted by their score:")
    print(sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names), 
             reverse=True))

def template_ens(model, X_dense, Y, names):
    print("running model", model)

    
    model.fit(X_dense, Y)
    print("Features sorted by their score:")
    ordering = sorted(zip(map(lambda x: round(x, 4), model.feature_importances_), names), 
             reverse=True)
    print(ordering) 

    rank = range(1, len(ordering)) 

    model_label = str(model).split("(")[0]



    print(rank)
    print(model_label) 


    df = pd.DataFrame(ordering, columns = ["score", "feature"])
    #rank = pd.DataFrame(range(0, len(ordering)), columns = ("rank"])
    df['model'] = model_label
    df.index.name="rank"
    df = df.reset_index()
    print(df)
    return df

def main():

    parser = argparse.ArgumentParser(description="Performs several different feature selections")
    parser.add_argument("--feature_names", action="store",
                                    help="Feature names in a file, format: a,b,c,d")
    parser.add_argument("--libsvm1_scale_file", action="store")
 
    parser.add_argument("--plainfeatmat", action="store")

    parser.add_argument("--output_file", action="store", required=True, default=None,
                                    help="Filename of output file")
    args = parser.parse_args()

    names = open(args.feature_names, "r").readline().rstrip("\n").split(",")
    print(names)
    if args.libsvm1_scale_file:
        X, Y = load_svmlight_file(args.libsvm1_scale_file)

        print(X)
        print(Y)

    

        #names = range(0, X.shape[1])

  
        print("change to dense")
        X_dense = X.toarray()


        print(X_dense)
   
    if args.plainfeatmat:
       
        df = pd.read_csv(args.plainfeatmat, sep = ",")

        X = df[names]
        X_dense = X.as_matrix()

        Y = df['label']

    final =pd.DataFrame()

    for model in (AdaBoostClassifier(), AdaBoostRegressor(), GradientBoostingClassifier(), GradientBoostingRegressor(), RandomForestClassifier(), RandomTreesEmbedding(), RandomForestRegressor(), ExtraTreesClassifier(), ExtraTreesRegressor()):

        final = final.append(template_ens(model, X_dense, Y, names))
        print(final)


    final.to_csv(args.output_file, index=False)
 

if __name__ == "__main__":
    main()





