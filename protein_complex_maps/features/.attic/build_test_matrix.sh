sbatch label_feature_matrix_test.sbatch 
sbatch select_labeled_1.sbatch 
sbatch select_labeled_0.sbatch 
sbatch convert_libsvm_test.sbatch 
cat libsvm_???_labeled_test_1 > libsvm_atobsc3_labeled_test_1
cat ???_labeled_test_1 > atobsc3_labeled_test_1
scp atobsc3_labeled_test_1 \claire@hopper.icmb.utexas.edu:/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/interaction_network
