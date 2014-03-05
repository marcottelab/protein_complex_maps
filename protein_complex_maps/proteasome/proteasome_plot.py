
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.read_data as rd

savefilename = "/home/kdrew/public_html/test/Proteasome_Hs_helaN_ph_hcw120_2.pdf"

sample_filename2 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome_sorted.txt"

proteosome_ids = ["ENSG00000128789","ENSG00000163636","ENSG00000013275","ENSG00000110801","ENSG00000041357","ENSG00000130706","ENSG00000108294","ENSG00000100567","ENSG00000126067","ENSG00000087191","ENSG00000106588","ENSG00000154611","ENSG00000099341","ENSG00000204264","ENSG00000095261","ENSG00000100764","ENSG00000165916","ENSG00000185627","ENSG00000100902","ENSG00000173692","ENSG00000100804","ENSG00000129084","ENSG00000100519","ENSG00000008018","ENSG00000101182","ENSG00000108671","ENSG00000108344","ENSG00000175166","ENSG00000143106","ENSG00000136930","ENSG00000132963","ENSG00000142507","ENSG00000197170","ENSG00000161057","ENSG00000159352","ENSG00000159377","ENSG00000115233","ENSG00000103035","ENSG00000101843"]

proteosome_lid_ids = ["ENSG00000099341", "ENSG00000108344", "ENSG00000197170", "ENSG00000108671", "ENSG00000163636", "ENSG00000185627"]

proteosome_core_ids = ["ENSG00000129084", "ENSG00000106588","ENSG00000100567", "ENSG00000041357", "ENSG00000143106", "ENSG00000100902","ENSG00000101182","ENSG00000154611","ENSG00000008018","ENSG00000100804","ENSG00000204264","ENSG00000126067","ENSG00000108294","ENSG00000159377","ENSG00000142507","ENSG00000136930"]

genename_map = {}
genename_map["ENSG00000008018"]="PSMB1"
genename_map["ENSG00000013275"]="PSMC4"
genename_map["ENSG00000041357"]="PSMA4"
genename_map["ENSG00000087191"]="PSMC5"
genename_map["ENSG00000095261"]="PSMD5"
genename_map["ENSG00000099341"]="PSMD8"
genename_map["ENSG00000100519"]="PSMC6"
genename_map["ENSG00000100567"]="PSMA3"
genename_map["ENSG00000100764"]="PSMC1"
genename_map["ENSG00000100804"]="PSMB5"
genename_map["ENSG00000100902"]="PSMA6"
genename_map["ENSG00000101182"]="PSMA7"
genename_map["ENSG00000101843"]="PSMD10"
genename_map["ENSG00000103035"]="PSMD7"
genename_map["ENSG00000106588"]="PSMA2"
genename_map["ENSG00000108294"]="PSMB3"
genename_map["ENSG00000108344"]="PSMD3"
genename_map["ENSG00000108671"]="PSMD11"
genename_map["ENSG00000110801"]="PSMD9"
genename_map["ENSG00000115233"]="PSMD14"
genename_map["ENSG00000126067"]="PSMB2"
genename_map["ENSG00000128789"]="PSMG2"
genename_map["ENSG00000129084"]="PSMA1"
genename_map["ENSG00000130706"]="ADRM1"
genename_map["ENSG00000132963"]="POMP"
genename_map["ENSG00000136930"]="PSMB7"
genename_map["ENSG00000142507"]="PSMB6"
genename_map["ENSG00000143106"]="PSMA5"
genename_map["ENSG00000154611"]="PSMA8"
genename_map["ENSG00000159352"]="PSMD4"
genename_map["ENSG00000159377"]="PSMB4"
genename_map["ENSG00000161057"]="PSMC2"
genename_map["ENSG00000163636"]="PSMD6"
genename_map["ENSG00000165916"]="PSMC3"
genename_map["ENSG00000173692"]="PSMD1"
genename_map["ENSG00000175166"]="PSMD2"
genename_map["ENSG00000185627"]="PSMD13"
genename_map["ENSG00000197170"]="PSMD12"
genename_map["ENSG00000204264"]="PSMB8"


#sample_file1 = open(sample_filename1, 'rb')
sample_file2 = open(sample_filename2, 'rb')

#data_matrix1, name_list1 = rd.read_datafile(sample_file1)
data_matrix2, name_list2 = rd.read_datafile(sample_file2)

#clean_data_matrix_pre_normalized1 = nu.remove_zero(data_matrix1)
clean_data_matrix_pre_normalized2 = nu.remove_zero(data_matrix2)

#clean_data_matrix_noised1 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized1)
clean_data_matrix_noised2 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized2)

#clean_data_matrix, name_list = rd.concat_data_matrix( clean_data_matrix_noised1, name_list1, clean_data_matrix_noised2, name_list2)
clean_data_matrix = clean_data_matrix_noised2
name_list = name_list2


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
		barcolor = "blue"
	if name_list[i] in proteosome_core_ids:
		barcolor = "red"

	data_subplots[i].bar(np.arange(len(data_row)), data_row, align='center', facecolor=barcolor)
	data_subplots[i].axes.set_yticklabels([],visible=False)
	data_subplots[i].set_ylabel("%s" % (genename_map[name_list[i]],),rotation='horizontal', color=barcolor, fontsize=8)

#kdrew: setting to show just a subset to highlight subcomplex
plt.xlim(40.0,80.0)
print plt.xlim()

plt.savefig(savefilename)


