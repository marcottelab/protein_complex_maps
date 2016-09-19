
#1. GroupID 2. Rank    3. Level 4. Species  5. ProteinID       6. evalue  7. QueryRange      8.ProteomeID     9. Hitlength        10. Sequence        11. Annotation

cd for_mySQL


#Need to delete geneid column, since not all have gene ids
#proteinID description, species
awk -F'\t' -OFS',' '{$5 $11 $4}' mySQL_orthology.tab > protein.csv
 
#proteinID groupID level
awk -F'\t' -OFS',' '{print $5 $1 $3}' mySQL_orthology.tab > protein_group.csv


#ProteinID Sequence
awk -F'\t' -OFS',' '{print $5 $10}' mySQL_orthology.tab > protein_sequence.csv


 
