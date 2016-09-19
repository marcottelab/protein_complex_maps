
BASEDIR=/project/kdrew/data/msblender_runs
=
#for exp in At_Col_0_indark_201505 wgSEC1  WheatGermSEC_07-2015_repr_unmapped_cdhit95 RiceL_IEX wgSEC1_repr_unmapped_cdhit95 WheatGermSEC_07-2015_uniprot wheatgermIEX WheatGermSEC_07-2015

for exp in Rice_201505_uniprot At_Col_0_indark_fraction_201504 At_Col_0_leaf_fraction_2014 At_Col_0_leaf_fraction_2015 
do
   rm /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv
   echo "ExperimentID,FractionID,Peptide,PeptideCount" >  /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv

   cd $BASEDIR/${exp}  
   for dir in *-Results
   do
       #echo "you are here"
       pwd
       cd $dir
       #echo $exp
       #echo $dir
       info=`echo $exp,$dir`
       echo $info
       tail -n +3 msblender.pep_count_mFDRpsm001 | awk -F'\t' -v OFS=',' -v info="$info" '{print info, $1, $2}'  >> /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv
       #rm msblender.pep_count_mFDRpsm001.tmp 
       cd ../
   done
   cd ../
done

