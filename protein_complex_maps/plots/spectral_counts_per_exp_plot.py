
import pickle as p
import glob as g
import numpy as np

import matplotlib as mlab
import matplotlib.pyplot as plt

tot_list = list()

#for filename in g.glob("/home/kdrew/data/protein_complex_maps/msds_pickles/fractions/*10.p"):
for filename in g.glob("/home/kdrew/data/protein_complex_maps/msds_pickles/Hs/fractions/*10.p"):
	f = open(filename,'rb')

	msds = p.load(f)

	dm = msds.get_data_matrix()

	tot_list.append(np.sum(dm))

	f.close()

#plt.hist(tot_list)
plt.bar(range(len(tot_list)), tot_list)

plt.show()


