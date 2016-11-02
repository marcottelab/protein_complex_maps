

BASEDIR=/MS.processed/Fusion_data/
for exp in OP_SelaginellaSEC_20160309

do
    consolidated=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/elutions/${exp}_elution.csv

    rm -f $consolidated    
    echo "ExperimentID,FractionID,Peptide,PeptideCount" >  $consolidated

    cd $BASEDIR/${exp}/ms1results  
    for frac in *.ms2.pep.csv
    do 
        tail -n +8 $frac > ${frac%.csv}.tmp
        inputfile=${frac%.csv}.tmp
 
        #If no peptides found in fraction, don't want to load into pandas 
        if ! [ -s "$inputfile" ] 
        then
            echo "no peptides from " $frac
            rm $inputfile
            continue
        fi
 


        echo "inputfile" $inputfile


        python /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/scripts/count_fusion_peptides.py $inputfile outfile.tmp

        rm $inputfile
        echo $exp
        fracID=${frac%.ms2.pep.csv}
        echo $fracID 
        info=`echo $exp,$fracID`
        echo $info
        awk -F'\t' -v OFS=',' -v info="$info" '{print info, $1, $2}' outfile.tmp  > outfile_info

        head outfile_info

        cat  $consolidated outfile_info > ${consolidated}.tmp
        echo ${consolidated}.tmp

        mv ${consolidated}.tmp $consolidated
        rm outfile.tmp
        rm outfile_info

   done
   cd ../
done
