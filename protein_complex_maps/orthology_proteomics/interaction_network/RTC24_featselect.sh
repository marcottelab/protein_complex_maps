INFILE=atobsc3_full_libsvm1.scale.txt
CONDITION=RTC24

awk -F' ' '{print $1, $22, $34, $36, $38, $40, $69, $70, $71, $73, $74, $79, $83, $85, $88, $92, $100, $101, $102, $104, $106, $108, $109, $113, $114}' $INFILE > ${INFILE%.txt}$CONDITION






