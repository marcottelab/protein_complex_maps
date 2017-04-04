INFILE=atobsc3_full_libsvm0.scaleByTrain.test
CONDITION=test_poissonplus

awk -F' ' '{print $1, $4, $8, $12, $16, $20, $27, $28, $33, $36, $40, $44, $46, $47, $50, $55, $58, $62, $66, $71, $74, $78, $84, $85, $92, $93, $98, $102, $106, $110, $114}' $INFILE > ${INFILE%.txt}_$CONDITION






