INFILE=atobsc3_full_libsvm1.scale.txt
CONDITION=arathtraesbraolselml

awk -F' ' '{print $1, $8, $12, $16, $20, $28, $36, $40, $44, $62, $66, $78, $92, $93}' $INFILE > ${INFILE%.txt}$CONDITION






