INFILE=atobsc3_full_libsvm0.scaleByTrain.test
CONDITION=test_spearman

awk -F' ' '{print $1, $3, $7, $11, $15, $19, $24, $25, $31, $35, $39, $43, $49, $54, $56, $60, $64, $70, $73, $77, $82, $83, $90, $91, $97, $101, $105, $109, $113}' $INFILE > ${INFILE%.txt}_$CONDITION






