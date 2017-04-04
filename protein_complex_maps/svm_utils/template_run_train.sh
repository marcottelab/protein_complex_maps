
PROJECT_DIR=/scratch/03491/cmcwhite/protein_complex_maps/protein_complex_maps
#WORKING_DIR=$PROJECT_DIR/orthology_proteomics/interaction_network
WORKING_DIR=$PROJECT_DIR/features
PROGRAM_DIR=/home1/03491/cmcwhite
BASENAME=$1
C=$2
G=$3
H=$4
PARAM=$5


#claire: train
#$PROGRAM_DIR/libsvm/svm-train -h $H -b 1  -c $C -g $G $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.txt $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.model_$PARAM

#/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/interaction_network/atobs_euNOG_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.txt

#claire: predict unlabeled set w/ test set on train model
#THIS now done in batch. Takes too long (multiple days)
#$PROGRAM_DIR/libsvm/svm-predict -b 1 $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain.txt $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.model_${PARAM}  $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob

#Prev already run

#################################### 0
#claire: probability ordered list of pairs for unlabeled and test set
#Is this ever used? no? Commented out for now
#python -u  /home/kdrew/scripts/protein_complex_maps/protein_complex_maps/features/svm_results2pairs.py --input_feature_matrix $PROJECT_DIR/interaction_network/atobs_euNOG_corumtrain_labeled.txt --input_results $PROJECT_DIR/interaction_network/atobs_euNOG_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob --output_file $PROJECT_DIR/interaction_network/atobs_euNOG_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_pairs_noself_nodups.txt --sep , --id_columns ID1 ID2

#claire: adding additional commands, previous commands ran successfully 6/12/16

#Get annotation columns from feature matrix
#THESE COLUMN NUMBERS LIKELY WRONG
#Takes about 5 minutes


#For OLD frozenset ID scheme
#awk -F',' '{print $(NF-8), $(NF-7), $(NF-2)}' ${BASENAME}_corumtrain_labeled.txt | tail -n +2 > tmp1

#get Label 0 IDs, then remove label column
#<1 minute
#For OLD frozenset ID scheme
#grep -e ' 0$' tmp1 | awk -F' ' '{print $1, $2}' > tmp1.0

#NEW PROCESS better


#This includes header
echo "start grep 1"
#grep -ve ',0.0$' ${BASENAME}_corumtrain_labeled.txt > ${BASENAME}_corumtrain_labeled1.txt

#This does not include header
echo "start grep 0"
#grep -e ',0.0$' ${BASENAME}_corumtrain_labeled.txt > ${BASENAME}_corumtrain_labeled0.txt


echo "cut out labels 1"
#awk -F',' '{print $1}' ${BASENAME}_corumtrain_labeled1.txt > ${BASENAME}_corumtrain_labels1.txt


#Both files have headers
echo "paste together 1"
#paste -d' ' ${BASENAME}_corumtrain_labels1.txt ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob > ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_labels.txt


echo "remove header 1"
#tail -n +2 ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_labels.txt >${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_labels_noheader.txt 



echo "get pos prob columns and label 1"
#awk '{print $1, $2, $4}' ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_labels.txt > ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_wprob.txt

echo "sort by column 3, remove duplicates"
#sort -n -k3 -r ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_wprob.txt > ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt


#Did grep -v "labels" atobsc_euNOG_corumtrain_labeled.libsvm0.scaleByTrain_c8_g05_h0.resultsWprob.cat > atobsc_euNOG_corumtrain_labeled.libsvm0.scaleByTrain_c8_g05_h0.resultsWprob

#This has no header
echo "cut out labels 0"
#awk -F',' '{print $1}' ${BASENAME}_corumtrain_labeled0.txt > ${BASENAME}_corumtrain_labels0.txt


#Remove header from resultsWprob
#echo "remove header 0"
#tail -n +2 ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob >${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_noheader 


echo "paste together 0"
paste -d' ' ${BASENAME}_corumtrain_labels0.txt ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob > ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_labels.txt

echo "get pos prob columns and label 0"
awk '{print $1, $2, $4}' ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_labels.txt > ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain.${PARAM}.resultsWprob_pairs_noself_wprob.txt

echo "sort by column 3, remove duplicates"
sort -n -k3 -r ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain.${PARAM}.resultsWprob_pairs_noself_wprob.txt > ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain.${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt



cat $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt  $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt | sort -g -k 3 -r > $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.libsvm0.scaleByTrain_${PARAM}.resultsWprob_pairs_noself_nodups_wprob_combined.txt

echo "done"









#~30 seconds
#For OLD frozenset ID scheme
#tail -n +2 ${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob  | awk -F' ' '{print $2}' > tmp2.0

#Join IDs with prob
#~15 seconds
#paste -d' ' tmp1.0 tmp2.0 > tmp3.0

#Takes 5-10 minutes
#sort -n -k3 -r tmp3.0 > tmp4.0

#mv tmp4.0 $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm0.scaleByTrain_${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt



##################### 1/-1
#claire: predict training set on train model 
#Already run
#$PROGRAM_DIR/libsvm/svm-predict -b 1 $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.txt $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.model_${PARAM}  $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob 




#claire: probability ordered list of pairs for training set with probability
#This doesn't work
#python $PROJECT_DIR/features/svm_results2pairs.py --input_feature_matrix $WORKING_DIR/${BASENAME}_corumtrain_labeled.txt --input_results $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob --output_file $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt --sep , --id_column ID --add_prob --label_not0 &

echo "done"
#get Labels -1/1 IDs, then remove label column
#20 seconds
#grep -ve ' 0$' tmp1 | awk -F' ' '{print $1, $2}' > tmp1.1

#instantaneous
#tail -n +2 ${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob  | awk -F' ' '{print $2}' > tmp2.1

#fast
#paste -d' ' tmp1.1 tmp2.1 > tmp3.1

#sort -n -k3 -r tmp3.1 > tmp4.1

#mv tmp4.1 $WORKING_DIR/${BASENAME}_corumtrain_labeled.libsvm1.scale.${PARAM}.resultsWprob_pairs_noself_nodups_wprob.txt


#claire: combine results from both test and training predictions (we could train everything at once so there is not this extra combination step but keeping consistent with previous work flow

#Takes like 15 minutes

