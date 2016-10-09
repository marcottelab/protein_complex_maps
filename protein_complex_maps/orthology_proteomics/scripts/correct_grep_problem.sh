for f in *virNOG*tophit.txt
do
   grep -P '\t1\tvirNOG' $f > ${f}.tmp
   mv ${f}.tmp $f

done














