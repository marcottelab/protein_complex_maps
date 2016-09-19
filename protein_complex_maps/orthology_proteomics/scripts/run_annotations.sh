BASEDIR=$( pwd )




#python $BASEDIR/scripts/annotate_tables.py annotation_files/KOG_annotations.csv correlation_elutions/atobs_euNOG_concat.csv identified_elutions/atobs/atobs_euNOG_labelled_concat.csv GroupID
#python $BASEDIR/scripts/annotate_tables.py annotation_files/all_annotations.csv correlation_elutions/atobs_euNOG_concat.csv identified_elutions/atobs/atobs_euNOG_labelled_concat.csv GroupID

python $BASEDIR/scripts/annotate_tables.py eggnog_output/arath.euNOG_map_annotations_simplified.tab identified_elutions/atobs/atobs_euNOG_labelled_concat.csv identified_elutions/atobs/atobs_euNOG_labelledarath_concat.csv GroupID


#python $BASEDIR/scripts/annotate_groups_w_arath.py eggnog_output/arath.euNOG_orthology.tab correlation_elutions/atobs_euNOG_concat.csv identified_elutions/atobs/atobs_euNOG_labelled_concat.csv


