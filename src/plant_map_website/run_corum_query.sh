
#Training set
#grep " .* .*"  allComplexesCore_geneid_merged06_testSplit1_humanOnly_under30.txt > corum_human_under_30_over_2.txt
CORUMLIST=$1
OUTFILE=$2
while read p
do

echo $p
#python protein_complex_maps/db_data/query_interactions_plmap.py --bait_complex $p --format gene_id --save_stats_file $OUTFILE --logname bla.txt
python get_distributions.py --bait_complex "$p" --format gene_id --save_stats_file $OUTFILE --logname bla.txt

done < $CORUMLIST


#Test set
#grep " .* .*" allComplexesCore_geneid_merged06_trainSplit_noTestOverlap_humanOnly.txt > /project/cmcwhite_protein_complex_maps/protein_complex_maps/db_data/corum_human_over2_test.txt

#while read p
#do

#echo $p >> corum_scores_text.txt
#python protein_complex_maps/db_data/query_interactions.py --bait_complex $p --format gene_id --save_stats_file corum_scores_test.txt --logname bla.txt

#done < protein_complex_maps/db_data/corum_human_over2_test.txt


