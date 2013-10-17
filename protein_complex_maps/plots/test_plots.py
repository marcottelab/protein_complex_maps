
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr


sample_filename = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"
#sample_filename = "/Users/kdrew/data/bborgeson/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"
#sample_filename = "/Users/kdrew/data/bborgeson/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome_exosc_randos.txt"
#sample_filename = "/Users/kdrew/data/bborgeson/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome.txt"

proteosome_ids = ["ENSG00000128789","ENSG00000163636","ENSG00000013275","ENSG00000110801","ENSG00000041357","ENSG00000130706","ENSG00000108294","ENSG00000100567","ENSG00000126067","ENSG00000087191","ENSG00000106588","ENSG00000154611","ENSG00000099341","ENSG00000204264","ENSG00000095261","ENSG00000100764","ENSG00000165916","ENSG00000185627","ENSG00000100902","ENSG00000173692","ENSG00000100804","ENSG00000129084","ENSG00000100519","ENSG00000008018","ENSG00000101182","ENSG00000108671","ENSG00000108344","ENSG00000175166","ENSG00000143106","ENSG00000136930","ENSG00000132963","ENSG00000142507","ENSG00000197170","ENSG00000161057","ENSG00000159352","ENSG00000159377","ENSG00000115233","ENSG00000103035","ENSG00000101843"]

exosome_ids = ["ENSG00000135698","ENSG00000178896","ENSG00000130713","ENSG00000123737","ENSG00000083520","ENSG00000120699","ENSG00000075914","ENSG00000223496","ENSG00000107371","ENSG00000171311","ENSG00000171824","ENSG00000077348"]

# example data
#mu = 100 # mean of distribution
#sigma = 15 # standard deviation of distribution
#x = mu + sigma * np.random.randn(10000)

sample_file = open(sample_filename,'rb')
#kdrew: eat header
line = sample_file.readline()
data = []
#for line in sample_file.readlines():
#	#print line
#	line_data = line.split()
#	#line_array = map(float,line_data)
#	print line_data
#	data.append(line_data)

name_list = []

for line in sample_file.readlines():
	#print line
	line_data = line.split()
	line_array = map(float,line_data[2:])
	name_list.append(line_data[0])
	print line_array
	data.append(line_array)
							

print len(data)


data_subplots = []
f, data_subplots = plt.subplots(len(data),1,sharex='col')

max_value = np.max(data)
print max_value

for i, data_row in enumerate(data):
	barcolor = "green"
	if name_list[i] in proteosome_ids:
		barcolor = "red"
	if name_list[i] in exosome_ids:
		barcolor = "blue"

	#data_subplots[i].bar(np.arange(len(data_row[2:])), map(float,data_row[2:]), align='center', facecolor=barcolor, alpha=0.5 )
	data_subplots[i].bar(np.arange(len(data_row)), data_row, align='center', facecolor=barcolor, alpha=0.5 )
	data_subplots[i].axes.set_yticklabels([],visible=False)
	#data_subplots[i].set_ylabel("%s:%s" % (i, data_row[0]),rotation='horizontal', color=barcolor)
	data_subplots[i].set_ylabel("%s" % (i,),rotation='horizontal', color=barcolor)
	data_subplots[i].axes.set_ylim(0,max_value)


corr_dict = dict()
corr_list = []
for data_row in data:
	for data_row2 in data:
		pr = pearsonr(map(float,data_row[2:]), map(float,data_row2[2:]))
		corr_list.append(pr[0])
		corr_dict[(data_row[0],data_row2[0])] = pr

print corr_dict

#pr = pearsonr(e_kd, r_be)
#print pr

# the histogram of the data
#n, bins, patches = plt.hist(line_array, num_bins, normed=0, facecolor='green', alpha=0.5)
#plt.bar(np.arange(len(line_array)), line_array, align='center', facecolor='green', alpha=0.5)
# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
#plt.xlabel('fraction')
#plt.ylabel('counts')
#plt.title(r'fraction')

# Tweak spacing to prevent clipping of ylabel
#plt.subplots_adjust(left=0.15)
plt.show()

num_bins = 100
n, bins, patches = plt.hist(corr_list, num_bins, normed=0, facecolor='green', alpha=0.5)
plt.show()

