INFILE=infile
CONDITION=condition

awk -F' ' '{print colselect}' $INFILE > ${INFILE%.txt}_$CONDITION






