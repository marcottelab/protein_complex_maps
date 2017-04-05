BASEDIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics


#cd $BASEDIR/identified_elutions/traes

#python $BASEDIR/scripts/combine_elution_protein.py traes proteins wgSEC1_raw_wide_elution_traes_proteins.csv wheatgermIEX_raw_wide_elution_traes_proteins.csv WheatGermSEC_07-2015_raw_wide_elution_traes_proteins.csv wgIEF_raw_wide_elution_traes_proteins.csv


cd $BASEDIR/annotated_elutions


#python $BASEDIR/scripts/annotate_tables.py $BASEDIR/annotation_files/uniprot_traes_protein_annotations.csv $BASEDIR/identified_elutions/traes/traes_proteins_concat.csv uniprot_annotated_traes_proteins_concat.csv Entry '\t' ',' uniprot



#python $BASEDIR/scripts/annotate_tables.py $BASEDIR/annotation_files/all_annotations.csv $BASEDIR/correlation_elutions/plants_euNOG_concat.csv euNOG_annotated_plant_euNOG_concat.csv GroupID ',' ',' euNOG_annot

for SPEC in arath orysj braol selml chlre
do

python $BASEDIR/scripts/annotate_tables.py $BASEDIR/eggnog_output/${SPEC}.euNOG_map_annotations.tab  $BASEDIR/annotated_elutions/euNOG_annotated_plant_euNOG_concat.csv species_${SPEC}_annotated_plant_euNOG_concat.csv GroupID '\t' ',' ${SPEC} --target_col AllMembers

done



python $BASEDIR/scripts/combine_elution_group.py all_species euNOG  $BASEDIR/eggnog_output/*euNOG_map_annotations.tab













