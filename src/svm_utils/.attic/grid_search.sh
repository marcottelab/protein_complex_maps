PROJECT_DIR=/project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/

input_libsvm1_scaled=$1
feature_columns=$2
feature_names=$3
#check it if exists


python ../../features/select_feature_columns.py --libsvm1_scaled $input_libsvm1_scaled --feature_columns $feature_columns --feature_header $feature_names



##Move this to the next step after feature selection
##claire: parameter sweep using training set
#features_RFC_arathtraesorysj_euNOG_corum_train_feature_selection.csv

echo "starting grid.py"
parallel -j8  python -u /home/kdrew/programs/libsvm-3.20/tools/grid.py --libsvm1_scaled {} ::: features_*${input_libsvm1_scaled}

#echo "parameter sweep by training set"
#python /home/kdrew/programs/libsvm-3.20/tools/grid.py ${EXP_ID}_corumtrain_labeled.libsvm1.scale.txt


echo "done"




