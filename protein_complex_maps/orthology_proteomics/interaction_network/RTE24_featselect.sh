INFILE=atobsc3_full_libsvm1.scale.txt
CONDITION=RTE24

awk -F' ' '{print $1, $3, $22, $30, $31, $34, $35, $42, $43, $46, $47, $50, $53, $62, $69, $74, $79, $80, $83, $92, $96, $98, $104, $108, $112}' $INFILE > ${INFILE%.txt}$CONDITION






