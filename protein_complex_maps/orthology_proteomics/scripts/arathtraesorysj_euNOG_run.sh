

PROJECT_DIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/

#kdrew: add label train
#python ~/scripts/protein_complex_maps/protein_complex_maps/features/add_label.py --input_feature_matrix $PROJECT_DIR/correlation_elutions/arathtraesorysj_euNOG_feature_matrix.txt --input_positives $PROJECT_DIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.train_ppis.txt --input_negatives $PROJECT_DIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.neg_train_ppis.txt --id_columns ID1 ID2 --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --fillna 0.0 --sep ,

#kdrew: add label test
#python ~/scripts/protein_complex_maps/protein_complex_maps/features/add_label.py --input_feature_matrix $PROJECT_DIR/correlation_elutions/arathtraesorysj_euNOG_feature_matrix.txt --input_positives  $PROJECT_DIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.test_ppis.txt --input_negatives  $PROJECT_DIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.neg_test_ppis.txt --id_columns ID1 ID2 --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_test_labeled.txt --fillna 0.0 --sep ,

#kdrew: failed after this point, rerunning just remaining commands

#kdrew: convert to libsvm format training set
python ~/scripts/protein_complex_maps/protein_complex_maps/features/feature2libsvm.py --input_feature_matrix $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.txt --features arath_euNOG_concat.txt.corr_poisson arathtraesorysj_euNOG_concat.txt.corr_poisson At_Col_0_indark_201505_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_indark_fraction_201504_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_leaf_fraction_2014_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_leaf_fraction_2015_raw_wide_elution_arath_euNOG.txt.corr_poisson orsyj_RiceL_IEX_raw_wide_elution_orysj_euNOG.csv_concat.txt.corr_poisson Rice_201505_uniprot_raw_wide_elution_orysj_euNOG.txt.corr_poisson RiceL_IEX_raw_wide_elution_orysj_euNOG.txt.corr_poisson traes_euNOG_concat.txt.corr_poisson wgSEC1_raw_wide_elution_traes_euNOG.txt.corr_poisson wheatgermIEX_raw_wide_elution_traes_euNOG.txt.corr_poisson WheatGermSEC_07-2015_raw_wide_elution_traes_euNOG.txt.corr_poisson  --label_column label --keep_labels -1 1 --sep ,

#kdrew: convert to libsvm format unlabeled set and test set
python ~/scripts/protein_complex_maps/protein_complex_maps/features/feature2libsvm.py --input_feature_matrix $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.txt --features arath_euNOG_concat.txt.corr_poisson arathtraesorysj_euNOG_concat.txt.corr_poisson At_Col_0_indark_201505_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_indark_fraction_201504_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_leaf_fraction_2014_raw_wide_elution_arath_euNOG.txt.corr_poisson At_Col_0_leaf_fraction_2015_raw_wide_elution_arath_euNOG.txt.corr_poisson orsyj_RiceL_IEX_raw_wide_elution_orysj_euNOG.csv_concat.txt.corr_poisson Rice_201505_uniprot_raw_wide_elution_orysj_euNOG.txt.corr_poisson RiceL_IEX_raw_wide_elution_orysj_euNOG.txt.corr_poisson traes_euNOG_concat.txt.corr_poisson wgSEC1_raw_wide_elution_traes_euNOG.txt.corr_poisson wheatgermIEX_raw_wide_elution_traes_euNOG.txt.corr_poisson WheatGermSEC_07-2015_raw_wide_elution_traes_euNOG.txt.corr_poisson  --label_column label --keep_labels 0 --sep ,

#kdrew: scale training
~/programs/libsvm-3.20/svm-scale -s $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale_parameters $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.txt > $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.txt 

#kdrew: scale unlabeled (w/ test set) by train’s scale parameters
~/programs/libsvm-3.20/svm-scale -r $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale_parameters $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.txt > $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.txt

##kdrew: parameter sweep using training set
python ~/programs/libsvm-3.20/tools/grid.py $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.txt 

