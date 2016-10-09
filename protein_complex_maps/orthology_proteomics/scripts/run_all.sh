spec=$1
level=$2
experimentID=$3
peptidepath=$4
BASEDIR=$( pwd )
contam=$5
raw_outfile_loc=$BASEDIR/correlation_elutions/$spec/

mkdir $raw_outfile_loc
#echo $BASEDIR

echo $spec
echo $level
echo $experimentID
echo $peptidepath



#Run once when there is a new experiment

#Pull experiment elutions from multiple files
#bash consolidate_MSblender_peptide_assignments.sh
#bash consolidate_MSblender_peptide_assignments_alt.sh

#Break all proteomes into peptides
#for f in proteomes/?????/*.fasta
#do

#   echo $f   
#   python scripts/trypsin.py --input $f --peptide_assignments ${f%.*}_peptides.csv --miss 2

#done


#This step gets tophit and nonhit eunog hits into one file to make _orthology.tab file
#bash $BASEDIR/scripts/get_eggnog_output2.sh $spec $level

#Place to put the peptide_assignments
mkdir $BASEDIR/peptide_assignments/${spec}

#mkdir identified_elutions/${spec}

#Lookup experimentally found peptides using orthogroups and individual proteins
python $BASEDIR/scripts/get_elution_ids.py $spec $level $experimentID $BASEDIR/eggnog_output/${spec}.${level}_orthology.tab $BASEDIR/elutions/${experimentID}_elution.csv $peptidepath $contam

python scripts/get_wideform_group.py identified_elutions/${spec}/${experimentID}_elution_${spec}_${level}.csv eggnog_output/${spec}.${level}_orthology.tab

python scripts/get_wideform_prot.py identified_elutions/${spec}/${experimentID}_elution_${spec}_proteins.csv

mv *.log logs/ 

cp identified_elutions/${spec}/${experimentID}_raw_wide*${level}.csv correlation_elutions


#python scripts/annotate_wideform.py identified_elutions/${spec}/${experimentID}_elution_${spec}_${level}.csv eggnog_output/${spec}.${level}_orthology.tab identified_elutions/${spec}/${experimentID}_alt_wide_elution_${spec}_proteins.csv $raw_outfile_loc






