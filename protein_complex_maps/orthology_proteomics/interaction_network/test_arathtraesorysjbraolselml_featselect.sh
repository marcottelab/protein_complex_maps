INFILE=xht_labeled_test_1
CONDITION=test_arathtraesorysjbraolselml

awk -F' ' '{print $1, $8, $12, $16, $20, $28, $36, $40, $44, $62, $66, $78, $92, $93}' $INFILE > ${INFILE%.txt}_$CONDITION






