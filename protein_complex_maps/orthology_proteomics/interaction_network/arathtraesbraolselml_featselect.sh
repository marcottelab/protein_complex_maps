INFILE=atobsc3_full_libsvm1.scale.oldmatched
CONDITION=arathtraesbraolselml_oldmatched_plants

awk -F' ' '{print $1, $8, $12, $16, $20, $28, $36, $40, $44, $58, $62, $66, $78, $92, $93, $114}' $INFILE > ${INFILE%.txt}$CONDITION






