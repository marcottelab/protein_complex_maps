BASEDIR=$( pwd )

spec=$1
level=$2
experiment=$3

cd $BASEDIR/peptide_assignments/${spec}/


#totpep_ortho=`wc -l nonredundant_orthogroup_peptides_${spec}_${level}.csv`
#totpep_prot=`wc -l nonredundant_protein_peptides_${spec}.csv`

id_ortho=`wc -l orthogroup_unique_peptides_${spec}_${level}.csv`
id_ortho=$(echo $id_ortho | cut -f1 -d' ')

id_prot=`wc -l protein_unique_peptides_${spec}.csv`
id_prot=$(echo $id_prot | cut -f1 -d' ')


cd $BASEDIR/identified_elutions/${spec}/

#identifications_ortho=`wc -l ${experiment}_elution_${spec}_${level}.csv`
#identifications_prot=`wc -l ${experiment}_elution_${spec}_proteins.csv`

spec_ortho=`awk -F',' '{print $5}' ${experiment}_elution_${spec}_${level}.csv| awk '{s+=$1} END {print s}'`
spec_prot=`awk -F',' '{print $5}' ${experiment}_elution_${spec}_proteins.csv| awk '{s+=$1} END {print s}'`


#echo $totpep_ortho
echo $id_ortho $level identifying_peptides $spec $experiment
echo $id_prot protein identifying_peptides $spec $experiment

#echo $identifications_ortho
echo $spec_ortho $level spectral_counts $spec $experiment



#echo $totpep_prot
#echo $identifications_prot
echo $spec_prot protein spectral_counts $spec $experiment

echo "- - - - -"
