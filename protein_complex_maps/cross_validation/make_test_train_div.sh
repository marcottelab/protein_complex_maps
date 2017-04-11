


awk -F',' '{print $1, $116}' atobsc3_labeled_full_1 > atobsc3_labels_labeled.txt 


echo "Separate positives and negatives so they can be selection in the same ratio"
grep " 1.0" atobsc3_labels_labeled.txt > atobsc3_labels_labeled_pos.txt
grep " -1.0" atobsc3_labels_labeled.txt > atobsc3_labels_labeled_neg.txt


echo "Divide sets into 5 cross validation sets (default)"
python train_test_divide.py --input_file atobsc3_labels_labeled_pos.txt 
python train_test_divide.py --input_file atobsc3_labels_labeled_neg.txt 

echo "Combine positives and negatives"
cat train_seg_1_atobsc3_labels_labeled_pos.txt train_seg_1_atobsc3_labels_labeled_neg.txt > train_seg_1_atobsc3_labels_labeled_posneg.txt
cat train_seg_2_atobsc3_labels_labeled_pos.txt train_seg_2_atobsc3_labels_labeled_neg.txt > train_seg_2_atobsc3_labels_labeled_posneg.txt
cat train_seg_3_atobsc3_labels_labeled_pos.txt train_seg_3_atobsc3_labels_labeled_neg.txt > train_seg_3_atobsc3_labels_labeled_posneg.txt
cat train_seg_4_atobsc3_labels_labeled_pos.txt train_seg_4_atobsc3_labels_labeled_neg.txt > train_seg_4_atobsc3_labels_labeled_posneg.txt
cat train_seg_5_atobsc3_labels_labeled_pos.txt train_seg_5_atobsc3_labels_labeled_neg.txt > train_seg_5_atobsc3_labels_labeled_posneg.txt

cat test_seg_1_atobsc3_labels_labeled_pos.txt test_seg_1_atobsc3_labels_labeled_neg.txt > test_seg_1_atobsc3_labels_labeled_posneg.txt
cat test_seg_2_atobsc3_labels_labeled_pos.txt test_seg_2_atobsc3_labels_labeled_neg.txt > test_seg_2_atobsc3_labels_labeled_posneg.txt
cat test_seg_3_atobsc3_labels_labeled_pos.txt test_seg_3_atobsc3_labels_labeled_neg.txt > test_seg_3_atobsc3_labels_labeled_posneg.txt
cat test_seg_4_atobsc3_labels_labeled_pos.txt test_seg_4_atobsc3_labels_labeled_neg.txt > test_seg_4_atobsc3_labels_labeled_posneg.txt
cat test_seg_5_atobsc3_labels_labeled_pos.txt test_seg_5_atobsc3_labels_labeled_neg.txt > test_seg_5_atobsc3_labels_labeled_posneg.txt



awk '{print $1, $2}' train_seg_1_atobsc3_labels_labeled_posneg.txt > train_seg_1_atobsc3_labels_posneg.txt
awk '{print $1, $2}' train_seg_2_atobsc3_labels_labeled_posneg.txt > train_seg_2_atobsc3_labels_posneg.txt
awk '{print $1, $2}' train_seg_3_atobsc3_labels_labeled_posneg.txt > train_seg_3_atobsc3_labels_posneg.txt
awk '{print $1, $2}' train_seg_4_atobsc3_labels_labeled_posneg.txt > train_seg_4_atobsc3_labels_posneg.txt
awk '{print $1, $2}' train_seg_5_atobsc3_labels_labeled_posneg.txt > train_seg_5_atobsc3_labels_posneg.txt


awk '{print $1, $2}' test_seg_1_atobsc3_labels_labeled_posneg.txt > test_seg_1_atobsc3_labels_posneg.txt
awk '{print $1, $2}' test_seg_2_atobsc3_labels_labeled_posneg.txt > test_seg_2_atobsc3_labels_posneg.txt
awk '{print $1, $2}' test_seg_3_atobsc3_labels_labeled_posneg.txt > test_seg_3_atobsc3_labels_posneg.txt
awk '{print $1, $2}' test_seg_4_atobsc3_labels_labeled_posneg.txt > test_seg_4_atobsc3_labels_posneg.txt
awk '{print $1, $2}' test_seg_5_atobsc3_labels_labeled_posneg.txt > test_seg_5_atobsc3_labels_posneg.txt














