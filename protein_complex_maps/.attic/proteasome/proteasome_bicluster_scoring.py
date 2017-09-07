import sys
import numpy as np
import protein_complex_maps.bicluster.bicluster as bc
import protein_complex_maps.normalization_util as nu
import protein_complex_maps.score_util as su
import protein_complex_maps.read_data as rd
import protein_complex_maps.random_sampling_util as rsu
import protein_complex_maps.plots.plot_bicluster as pb



np.random.seed(1234)
sample_filename1 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_hekN_1108_psome_exosc_randos.txt"
sample_filename2 = "/home/kdrew/data/protein_complex_maps/sample_data/Hs_helaN_ph_hcw120_2_psome_exosc_randos.txt"

sample_file1 = open(sample_filename1, 'rb')
sample_file2 = open(sample_filename2, 'rb')

data_matrix1, name_list1 = rd.read_datafile(sample_file1)
data_matrix2, name_list2 = rd.read_datafile(sample_file2)

#kdrew: remove columns with all zeros
clean_data_matrix_pre_normalized1 = nu.remove_zero(data_matrix1)
clean_data_matrix_pre_normalized2 = nu.remove_zero(data_matrix2)

clean_data_matrix_noised1 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized1)
clean_data_matrix_noised2 = nu.add_noise_over_columns(clean_data_matrix_pre_normalized2)

clean_data_matrix_normalized1 = nu.normalize_over_columns( clean_data_matrix_noised1 )
clean_data_matrix_normalized2 = nu.normalize_over_columns( clean_data_matrix_noised2 )

clean_data_matrix_min2one1 = nu.min_to_one_scale( clean_data_matrix_normalized1 )
clean_data_matrix_min2one2 = nu.min_to_one_scale( clean_data_matrix_normalized2 )

clean_data_matrix, name_list = rd.concat_data_matrix( clean_data_matrix_noised1, name_list1, clean_data_matrix_noised2, name_list2)
print clean_data_matrix

clean_data_matrix_normalized, name_list_normalized = rd.concat_data_matrix( clean_data_matrix_normalized1, name_list1, clean_data_matrix_normalized2, name_list2)

clean_data_matrix_min2one, name_list_min2one= rd.concat_data_matrix( clean_data_matrix_min2one1, name_list1, clean_data_matrix_min2one2, name_list2)

mdn_randsamp = rsu.RandomSampling( su.multiple_dot, sample_module = np.random ) 
mdpu_randsamp = rsu.RandomSampling( su.multiple_dot_per_unit, sample_module = np.random ) 
#sum_randsamp = rsu.RandomSampling( su.sum_cells, sample_module = np.random ) 

#kdrew: 10 (Rpn12) is part of the lid but has little signal in these columns, 32 (Rpn11) and 61 (Rpn8) are enclosed by lid
proteasome_lid_bicluster_compact = bc.Bicluster(rows = [25,34,37,40,48,32,61], cols = [41,123], random_module=np.random)
proteasome_lid_bicluster = bc.Bicluster(rows = [25,34,37,40,48,32,61], cols = [40,41,42,122,123], random_module=np.random)
proteasome_lid_bicluster_w_common_column = bc.Bicluster(rows = [25,34,37,40,48,32,61], cols = [40,41,42,122,123,133,134,135], random_module=np.random)
proteasome_lid_bicluster_background = bc.Bicluster(rows = [3,4,6,7,8,25,34,37,40,48,32,61], cols = [40,41,42,122,123,133,134,135], random_module=np.random)

proteasome_lid_bicluster_submatrix = bc.Bicluster(rows = [5,6,7,8,9,10,11], cols = [0,1,2,3,4], random_module=np.random)

#pb.plot_bicluster(clean_data_matrix, proteasome_lid_bicluster, savefilename="/home/kdrew/public_html/test/bicluster_proteasome_lid_plot.pdf" )
pb.plot_bicluster(proteasome_lid_bicluster_background.get_submatrix(clean_data_matrix), proteasome_lid_bicluster_submatrix, savefilename="/home/kdrew/public_html/test/bicluster_proteasome_lid_plot.pdf" )

print "\n"

print "clean_data_matrix: normal bicluster"
proteasome_lid_bicluster.print_submatrix(clean_data_matrix)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster )
#sum_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster )

print "\n"

print "clean_data_matrix: wide bicluster"
proteasome_lid_bicluster_w_common_column.print_submatrix(clean_data_matrix)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster_w_common_column )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster_w_common_column )
#sum_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster_w_common_column )

print "\n"

print "clean_data_matrix: compact bicluster"
proteasome_lid_bicluster_compact.print_submatrix(clean_data_matrix)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster_compact )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix, proteasome_lid_bicluster_compact )

print "\n"

print "normalized scaled data: normal bicluster"
proteasome_lid_bicluster.print_submatrix(clean_data_matrix_min2one)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster )

print "\n"

print "normalized scaled data: wide bicluster"
proteasome_lid_bicluster_w_common_column.print_submatrix(clean_data_matrix_min2one)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster_w_common_column )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster_w_common_column )

print "\n"

print "normalized scaled data: compact bicluster"
proteasome_lid_bicluster_compact.print_submatrix(clean_data_matrix_min2one)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster_compact)
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_min2one, proteasome_lid_bicluster_compact)

print "\n"

print "normalized data: normal bicluster"
proteasome_lid_bicluster.print_submatrix(clean_data_matrix_normalized)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster )

print "\n"

print "normalized data: wide bicluster"
proteasome_lid_bicluster_w_common_column.print_submatrix(clean_data_matrix_normalized)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster_w_common_column )
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster_w_common_column )

print "\n"

print "normalized data: compact bicluster"
proteasome_lid_bicluster_compact.print_submatrix(clean_data_matrix_normalized)
print "\tmultidot"
mdn_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster_compact)
print "\tmultidot perunit"
mdpu_randsamp.print_score( clean_data_matrix_normalized, proteasome_lid_bicluster_compact)
