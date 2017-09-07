
commands_file=atobsc3_full_libsvm1.scale.poissonplus.pr_COMMANDS.sh
rm $commands_file

for segnum in 1 2 3 4 5
do
    for cvalue in 0.00390625 0.0078125 2 32 128 250 500 1000 
    do
    
        trains=`ls atobsc3_full_libsvm1.scale.poissonplus.train.seg${segnum}.c${cvalue}*resultsWprob_pairs_sort | tr '\n' ' '` 
        tmp_labels_train=${trains//atobsc3_full_libsvm1.scale.poissonplus.train.}
        labels_train=${tmp_labels_train//.resultsWprob_pairs_sort}
        
        tests=`ls atobsc3_full_libsvm1.scale.poissonplus.test.seg${segnum}.c${cvalue}*resultsWprob_pairs_sort| tr '\n' ' '` 
        tmp_labels_test=${tests//atobsc3_full_libsvm1.scale.poissonplus.test.}
        labels_test=${tmp_labels_test//.resultsWprob_pairs_sort}
        
        
        echo "python /work/03491/cmcwhite/cross_validation/scripts/prcurve.py --results_wprob $tests --input_positives /work/03491/cmcwhite/cross_validation/corum/test_seg_${segnum}_atobsc3_labels_labeled_pos.labels --input_negatives /work/03491/cmcwhite/cross_validation/corum/test_seg_${segnum}_atobsc3_labels_labeled_neg.labels --output_file  atobsc3_full_libsvm1.scale.poissonplus.test.seg${segnum}.c${cvalue//.}_g${gvalue//.}.resultsWprob.pdf --labels $labels_test" >> $commands_file
        
        
        echo "python /work/03491/cmcwhite/cross_validation/scripts/prcurve.py --results_wprob $trains --input_positives /work/03491/cmcwhite/cross_validation/corum/train_seg_${segnum}_atobsc3_labels_labeled_pos.labels --input_negatives /work/03491/cmcwhite/cross_validation/corum/train_seg_${segnum}_atobsc3_labels_labeled_neg.labels --output_file  atobsc3_full_libsvm1.scale.poissonplus.train.seg${segnum}.c${cvalue//.}_g${gvalue//.}.resultsWprob.pdf --labels $labels_train" >> $commands_file
    done
done
