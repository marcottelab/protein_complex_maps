INFILE=atobsc3_full_libsvm0.scaleByTrain.test
CONDITION=test_arathtraesorysjbraolselml_plants


awk -F' ' '{print $1, $8, $12, $16, $20, $28, $36, $40, $44, $58, $62, $66, $78, $92, $93, $114}' $INFILE > ${INFILE%.txt}$CONDITION







