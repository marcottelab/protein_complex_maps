


#BASEDIR=/project/kdrew/data/msblender_runs
#for exp in At_Col_0_indark_201505 wgSEC1  WheatGermSEC_07-2015_repr_unmapped_cdhit95 RiceL_IEX wgSEC1_repr_unmapped_cdhit95 WheatGermSEC_07-2015_uniprot wheatgermIEX WheatGermSEC_07-2015

#for exp in Rice_201505_uniprot At_Col_0_indark_fraction_201504 At_Col_0_leaf_fraction_2014 At_Col_0_leaf_fraction_2015 
#BASEDIR=/MS/processed/Lumos
BASEDIR=/MS/processed/Fusion_data
#
   # for exp in OP_SelaginellaSEC_20160309 Fern_Frond_WWC_20151217-20160119
#for exp in Broccolinuclei_6-2016
#for exp in Fern_Frond_WWC_20151217-20160119
#for exp in selaginella_WWC wgIEF BroccoliNE_WWC
#for exp in wgIEF
#for exp in selaginella_WWC
#for exp in BroccoliNE_WWC selaginella_WWC
#for exp in Chlamydomonas_WWC_9-2016
#for exp in Broccolinuclei_IEF_9-2016
for exp in OP_Fernfrond_SEC_112016 
do
   #rm /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv
   echo "ExperimentID,FractionID,Peptide,PeptideCount" >  /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv

   cd $BASEDIR/${exp}/msblender
   for pep in *.pep_count_mFDRpsm001
   do

       #echo "you are here"
       #pwd
       echo $exp
   
       frac=${pep%.pep_count_mFDRpsm001}
       #echo $dir
       info=`echo $exp,$frac`
       echo $info
       tail -n +3 $pep | awk -F'\t' -v OFS=',' -v info="$info" '{print info, $1, $2}'  >> /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv
       #rm msblender.pep_count_mFDRpsm001.tmp 
   done
   cd ../
done

