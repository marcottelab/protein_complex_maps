#bash scripts/consolidate_MSblender_output.sh pep_counts_folder/Hs_hekN_1108 Hs_hekN_1108 elutions
#/pep_counts_folder/
pep_count_dir=$1
expID=$2
#/elutions
elution_dir=$3





echo "ExperimentID,FractionID,Peptide,PeptideCount" >  $elution_dir/${expID}_elution.csv

for pep in $pep_count_dir/*.pep_count_mFDRpsm001
do

    #echo "you are here"
    #pwd
    #echo $expID
    pep_name=${pep##*/}
  
    frac=${pep_name%.pep_count_mFDRpsm001}
    #echo $dir
    info=`echo $expID,$frac`
    echo $info
    tail -n +3 $pep | awk -F'\t' -v OFS=',' -v info="$info" '{print info, $1, $2}'  >> $elution_dir/${expID}_elution.csv
done

