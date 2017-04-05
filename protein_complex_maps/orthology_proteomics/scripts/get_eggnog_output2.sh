echo "combine the nonhits and tophit files" 


spec=$1
level=$2

echo "concatenating peptide lists"


#Fix for this will be making the output of eggnog be comma separated

cd eggnog_output

tail -n +2 ${spec}.${level}_tophit.txt | cat > ${spec}.${level}_orthology.tmp
tail -n +2 ${spec}.${level}_nonhits.txt | cat >> ${spec}.${level}_orthology.tmp

echo -e "GroupID\tRank\tLevel\tSpecies\tProteinID\tevalue\tQueryRange\tProteomeID\tHitlength\tSequence\tAnnotation" > ${spec}.${level}_orthology.tmp2


sort -u ${spec}.${level}_orthology.tmp >> ${spec}.${level}_orthology.tmp2


#Add new annotations to original pipeline on stampede, but for now...


python ../scripts/annotate_orthology_w_KOG.py ../annotation_files/KOG_annotations.tsv ${spec}.${level}_orthology.tmp2 ${spec}.${level}_orthology.tab






rm ${spec}.${level}_orthology.tm* 


