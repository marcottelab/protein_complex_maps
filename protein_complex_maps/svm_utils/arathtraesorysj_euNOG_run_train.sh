

PROJECT_DIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/

#kdrew: train (grid search not finished but reasonable parameters)
#~/programs/libsvm-3.20/svm-train -b 1  -c 32.0 -g 0.0078125 $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.model_c32_g0078125

#kdrew: predict unlabeled set w/ test set on train model
#~/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.txt $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.model_c32_g0078125  $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob

#kdrew: probability ordered list of pairs for unlabeled and test set
#python ~/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --input_results $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns ID1 ID2

#kdrew: adding additional commands, previous commands ran successfully 6/12/16

#kdrew: probability ordered list of pairs for unlabeled and test set with probability
python ~/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --input_results $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob

#kdrew: predict training set on train model
~/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.model_c32_g0078125  $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob

#kdrew: probability ordered list of pairs for training set with probability
python ~/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.txt --input_results $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob --output_file $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob --label_not0


#kdrew: combine results from both test and training predictions (we could train everything at once so there is not this extra combination step but keeping consistent with previous work flow
cat $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.resultsWprob_pairs_noself_nodups_wprob.txt  $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt |sort -g -k 3 -r > $PROJECT_DIR/interaction_network/arathtraesorysj_euNOG_corum_train_labeled.libsvm1.scale.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob_combined.txt


