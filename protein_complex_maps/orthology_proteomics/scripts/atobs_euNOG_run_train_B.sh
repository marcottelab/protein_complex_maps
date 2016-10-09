#From Grid.py log2c=1 log2g=1 rate=71.30..

PROJECT_DIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics

#claire: train
/home/kdrew/programs/libsvm-3.20/svm-train -h 0 -b 1  -c 2 -g 2 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g2_h0


#claire: predict unlabeled set w/ test set on train model
/home/kdrew/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g2_h0  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob

#claire: probability ordered list of pairs for unlabeled and test set
python -u  /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns ID1 ID2

#claire: adding additional commands, previous commands ran successfully 6/12/16

#claire: probability ordered list of pairs for unlabeled and test set with probability
python -u /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob

#claire: predict training set on train model 
/home/kdrew/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g2_h0  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g2_h0.resultsWprob


#claire: probability ordered list of pairs for training set with probability
python -u /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g2_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g2_h0.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob --label_not0


#claire: combine results from both test and training predictions (we could train everything at once so there is not this extra combination step but keeping consistent with previous work flow
cat $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g2_h0.resultsWprob_pairs_noself_nodups_wprob.txt  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob_pairs_noself_nodups_wprob.txt |sort -g -k 3 -r > $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.libsvm0.scaleByTrain_c2_g2_h0.resultsWprob_pairs_noself_nodups_wprob_combined.txt
