
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import protein_complex_maps.read_data as rd
import protein_complex_maps.normalization_util as nu


sample_filename1 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"
sample_filename2 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome_exosc_randos.txt"

proteosome_ids = ["ENSG00000128789","ENSG00000163636","ENSG00000013275","ENSG00000110801","ENSG00000041357","ENSG00000130706","ENSG00000108294","ENSG00000100567","ENSG00000126067","ENSG00000087191","ENSG00000106588","ENSG00000154611","ENSG00000099341","ENSG00000204264","ENSG00000095261","ENSG00000100764","ENSG00000165916","ENSG00000185627","ENSG00000100902","ENSG00000173692","ENSG00000100804","ENSG00000129084","ENSG00000100519","ENSG00000008018","ENSG00000101182","ENSG00000108671","ENSG00000108344","ENSG00000175166","ENSG00000143106","ENSG00000136930","ENSG00000132963","ENSG00000142507","ENSG00000197170","ENSG00000161057","ENSG00000159352","ENSG00000159377","ENSG00000115233","ENSG00000103035","ENSG00000101843"]

proteosome_lid_ids = ["ENSG00000099341", "ENSG00000108344", "ENSG00000197170", "ENSG00000108671", "ENSG00000163636", "ENSG00000185627"]

proteosome_core_ids = ["ENSG00000129084", "ENSG00000106588","ENSG00000100567", "ENSG00000041357", "ENSG00000143106", "ENSG00000100902","ENSG00000101182","ENSG00000154611","ENSG00000008018","ENSG00000100804","ENSG00000204264","ENSG00000126067","ENSG00000108294","ENSG00000159377","ENSG00000142507","ENSG00000136930"]

exosome_ids = ["ENSG00000135698","ENSG00000178896","ENSG00000130713","ENSG00000123737","ENSG00000083520","ENSG00000120699","ENSG00000075914","ENSG00000223496","ENSG00000107371","ENSG00000171311","ENSG00000171824","ENSG00000077348"]

sample_file1 = open(sample_filename1, 'rb')
sample_file2 = open(sample_filename2, 'rb')

data_matrix1, name_list1 = rd.read_datafile(sample_file1)
data_matrix2, name_list2 = rd.read_datafile(sample_file2)

clean_data_matrix_pre_normalized1 = nu.remove_zero(data_matrix1)
clean_data_matrix_pre_normalized2 = nu.remove_zero(data_matrix2)

clean_data_matrix_noised1 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized1)
clean_data_matrix_noised2 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized2)

clean_data_matrix, name_list = rd.concat_data_matrix( clean_data_matrix_noised1, name_list1, clean_data_matrix_noised2, name_list2)

##kdrew: eat header
#line = sample_file.readline()
#data = []
##for line in sample_file.readlines():
##	#print line
##	line_data = line.split()
##	#line_array = map(float,line_data)
##	print line_data
##	data.append(line_data)
#
#name_list = []
#
#for line in sample_file.readlines():
#	#print line
#	line_data = line.split()
#	line_array = map(float,line_data[2:])
#	name_list.append(line_data[0])
#	print line_array
#	data.append(line_array)
							

data_subplots = []
f, data_subplots = plt.subplots(len(clean_data_matrix),1,sharex='col')

max_value = np.max(clean_data_matrix)
print max_value

for i, data_r in enumerate(clean_data_matrix):
	print "pos: %s name: %s" % (i, name_list[i])
	data_row = np.array(data_r.reshape(-1))[0]
	barcolor = "gray"
	if name_list[i] in proteosome_ids:
		barcolor = "green"
	if name_list[i] in proteosome_lid_ids:
		barcolor = "purple"
	if name_list[i] in proteosome_core_ids:
		barcolor = "red"
	if name_list[i] in exosome_ids:
		barcolor = "blue"

	#data_subplots[i].bar(np.arange(len(data_row[2:])), map(float,data_row[2:]), align='center', facecolor=barcolor, alpha=0.5 )
	data_subplots[i].bar(np.arange(len(data_row)), data_row, align='center', facecolor=barcolor, alpha=0.5 )
	data_subplots[i].axes.set_yticklabels([],visible=False)
	#data_subplots[i].set_ylabel("%s:%s" % (i, data_row[0]),rotation='horizontal', color=barcolor)
	data_subplots[i].set_ylabel("%s" % (i,),rotation='horizontal', color=barcolor)
	#data_subplots[i].axes.set_ylim(0,max_value)


plt.show()


#corr_dict = dict()
#corr_list = []
#for data_row in data:
#	for data_row2 in data:
#		pr = pearsonr(map(float,data_row[2:]), map(float,data_row2[2:]))
#		corr_list.append(pr[0])
#		corr_dict[(data_row[0],data_row2[0])] = pr
#
#print corr_dict

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

#num_bins = 100
#n, bins, patches = plt.hist(corr_list, num_bins, normed=0, facecolor='green', alpha=0.5)
#plt.show()

