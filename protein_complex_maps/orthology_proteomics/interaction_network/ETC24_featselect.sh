INFILE=atobsc3_full_libsvm1.scale.txt
CONDITION=ETC24

awk -F' ' '{print $1, $34, $36, $38, $40, $47, $48, $68, $69, $73, $74, $78, $79, $82, $83, $84, $92, $102, $106, $108, $109, $110, $112, $113, $114}' $INFILE > ${INFILE%.txt}$CONDITION






