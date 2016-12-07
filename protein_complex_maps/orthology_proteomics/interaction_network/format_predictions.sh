

# format test

awk -F',' '{print $1}' atobsc3_labeled_test_1 > atobsc3_labeled_test_1_labels


#/home/claire/Programs/libsvm/svm-scale -r atobsc3_full.libsvm1.scale_parameters libsvm_atobsc3_labeled_test_1 > atobsc3_full_test_libsvm1.scaleByTrain.txt


#/home/claire/Programs/libsvm/svm-predict -b 1 atobsc3_full_test_libsvm1.scaleByTrain.txt atobsc3_full_libsvm1.scale.16thread.scale.model_h0_c32_g05 atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob

paste -d' ' atobsc3_labeled_test_1_labels atobsc3_full_test_scaleByTrain_h0_c32_g05.resultsWprob > atobsc3_full_test_scaleByTrain_h0_c32_g05.resultsWprob_labels

#tail -n +2 atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_labels > atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_labels_noheader

#awk -F' ' '{print $1, $2, $4}' atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_labels_noheader > atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_pairs

#sort -g -k3 -r atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_pairs > atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_pairs_sort


#nohup python ../../features/prcurve.py --results_wprob old/arathtraesorysjbraolselml_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt atobsc3_full_test_libsvm1.scaleByTrain_h0_c32_g05.resultsWprob_pairs_sort --input_positives ../corum/ordered_nonredundant_allComplexesCore_mammals_euNOG_merged06.test_ppis.txt --input_negatives ../corum/ordered_nonredundant_allComplexesCore_mammals_euNOG_merged06.neg_test_ppis.txt --output_file atobsc3_c32_g05_n.prcurve.pdf --labels arathtraesorysjbraolselml atobsc3 & 




#replace space with tab















