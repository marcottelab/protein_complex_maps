INFILE=atobsc3_full_libsvm1.scale.txt
CONDITION=ABC24

awk -F' ' '{print $1, $9, $13, $22, $24, $37, $38, $48, $53, $55, $77, $78, $79, $80, $82, $83, $84, $85, $87, $88, $92, $99, $102, $103, $112}' $INFILE > ${INFILE%.txt}$CONDITION






