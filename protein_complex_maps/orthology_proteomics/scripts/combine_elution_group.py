import sys
import pandas as pd
def join_elutions(elut1, elut2):

    elut1 = elut1.set_index(['GroupID'])
    elut2 = elut2.set_index(['GroupID'])
    final = elut1.join(elut2, how = "outer")
    final = final.reset_index()
    return final


firstname = sys.argv[3]
final = pd.read_csv(firstname, index_col=False)   

for f in sys.argv[4:]:
   f = pd.read_csv(f, index_col=False)
   final = join_elutions(final, f)

outfilename = sys.argv[1]+ "_" + sys.argv[2] + "_" + "concat.csv"

final = final.fillna(0)
final.to_csv(outfilename, index=False)


