#From Grid.py log2c=1 log2g=3 rate=71.9979

PROJECT_DIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics

#claire: train
/home/kdrew/programs/libsvm-3.20/svm-train -h 0 -b 1  -c 2 -g 8 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g8_h0

#/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.txt

#claire: predict unlabeled set w/ test set on train model
/home/kdrew/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g8_h0  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob


#################################### 0
#claire: probability ordered list of pairs for unlabeled and test set
#Is this ever used? no? Commented out for now
#python -u  /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns ID1 ID2

#claire: adding additional commands, previous commands ran successfully 6/12/16

#Get annotation columns from feature matrix
#Takes about 5 minutes
awk -F',' '{print $(NF-8), $(NF-7), $(NF-2)}' atobs_euNOG_corum_train_downsamp_labeled.txt | tail -n +2 > tmp1

#claire: probability ordered list of pairs for unlabeled and test set with probabilityi
#This step hasn't successfully run yet

#get Label 0 IDs, then remove label column
#<1 minute
grep -e ' 0$' tmp1 | awk -F' ' '{print $1, $2}' > tmp1.0

#~30 seconds
tail -n +2 atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob  | awk -F' ' '{print $2}' > tmp2.0

#Join IDs with prob
#~15 seconds
paste -d' ' tmp1.0 tmp2.0 > tmp3.0

#Takes 5-10 minutes
sort -n -k3 -r tmp3.0 > tmp4.0

mv tmp4.0 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt


#This takes too much memory to run on Hopper, replaced with the above Bash
#python -u /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob

##################### 1/-1
#claire: predict training set on train model 
/home/kdrew/programs/libsvm-3.20/svm-predict -b 1 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.txt $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.model_c2_g8_h0  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob 


#claire: probability ordered list of pairs for training set with probability
nohup python -u /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_columns ID1 ID2 --add_prob --label_not0 &

#get Labels -1/1 IDs, then remove label column
#20 seconds
grep -ve ' 0$' tmp1 | awk -F' ' '{print $1, $2}' > tmp1.1

#instantaneous
tail -n +2 atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob  | awk -F' ' '{print $2}' > tmp2.1

#fast
paste -d' ' tmp1.1 tmp2.1 > tmp3.1

sort -n -k3 -r tmp3.1 > tmp4.1

mv tmp4.1 $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt


#claire: combine results from both test and training predictions (we could train everything at once so there is not this extra combination step but keeping consistent with previous work flow

#Takes like 15 minutes

cat $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt  $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob.txt | sort -g -k 3 -r > $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale.libsvm0.scaleByTrain_c2_g8_h0.resultsWprob_pairs_noself_nodups_wprob_combined.txt
