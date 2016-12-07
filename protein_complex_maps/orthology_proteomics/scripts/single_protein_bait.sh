
#Need to rerun human to get nonhits
#python ../scripts/annotate_tables.py ../eggnog_output/human.euNOG_entry_tophit.txt Entry_jbw_dynein_baits.txt Entry_euNOG_jbw_dynein_baits.txt Entry '\t' ',' euNOG --target_col GroupID

#awk -F',' '{print $2}' Entry_euNOG_jbw_dynein_baits.txt > euNOG_jbw_dynein_baits.txt

FEATURE_MATRIX=chlre_euNOG_concat.txt.corr_poisson.pairs.ordered
echo $FEATURE_MATRIX
while read BAIT
do
#Use scaled feature matrix?
echo $BAIT
grep $BAIT $FEATURE_MATRIX > ${BAIT}.tmp

#Sort by each selected feature, take top 10
for feature_column in 5 6 7 8 9 
do

   sort -n -k$feature_column ${BAIT}.tmp > ${BAIT}.tmp2
   awk -F',' -OFS=',' '{print $1, $feature_column}' ${BAIT}.tmp2 > ${BAIT}.f${feature_column}.txt
   
done


done < euNOG_jbw_dynein_baits.txt











