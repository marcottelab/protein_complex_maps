import sys
import pandas as pd
def join_elutions(elut1, elut2):
    #print elut1
    #elut1 = elut1.set_index(['GroupID', 'ProteinID'])
    #print elut1
    #elut2 = elut2.set_index(['GroupID', 'ProteinID'])
    #print elut2
    elut1['key'] = elut1[['GroupID', 'ProteinID']].apply(lambda x: ' '.join(x), axis=1)
    elut2['key'] = elut2[['GroupID', 'ProteinID']].apply(lambda x: ' '.join(x), axis=1)
    elut1 = elut1.set_index(['key'])
    elut2 = elut2.set_index(['key'])
    final = elut1.merge(elut2, how = "outer")
    final = final.reset_index(drop=True)
    #final = final.reset_index()
    #final  = final.rename(columns={'level_0':'GroupID'})
    #final  = final.rename(columns={'index':'ProteinID'})
    return final


firstname = sys.argv[2]
final = pd.read_csv(firstname, index_col=False)   


for f in sys.argv[3:]:
   print f
   fopen = pd.read_csv(f, index_col=False,sep=",")
   
   final = join_elutions(final, fopen)
   print final.columns.values
outfilename = sys.argv[1]

final = final.fillna(0)
print final
print outfilename

final.to_csv(outfilename, index=False)



