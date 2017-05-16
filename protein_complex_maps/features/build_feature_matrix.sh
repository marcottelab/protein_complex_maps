


FEATURES=$1
EXP_ID=$2
EXISTING_MATRIX=$3
FEATURE_STR=$(cat $FEATURES | tr '\n' ' ')

echo ID: $EXP_ID
echo Feature file : $FEATURES
echo Features: $FEATURE_STR
echo Existing matrix: $EXISTING_MATRIX
#nohup python -u  /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files $FEATURE_STR --output_file ${EXP_ID}_feature_matrix.txt --prev_feature_matrix $EXISTING_MATRIX >  nohup_${EXP_ID}_feature_matrix.out&
python   /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files $FEATURE_STR --output_file ${EXP_ID}_feature_matrix.txt --prev_feature_matrix $EXISTING_MATRIX --store_interval 20 



if [ ! -f ${EXP_ID}_feature_matrix.txt ]; then
    echo "Building feature matrix $EXP_ID with features $FEATURE_STR failed :("
fi


