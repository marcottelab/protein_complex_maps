

EXP_ID=$1
MATRIX1=$2
MATRIX2=$3
echo ID: $EXP_ID
echo $MATRIX1
echo $MATRIX2

#nohup python -u  /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files $FEATURE_STR --output_file ${EXP_ID}_feature_matrix.txt --prev_feature_matrix $EXISTING_MATRIX >  nohup_${EXP_ID}_feature_matrix.out&
python /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/combine_feature_matrices.py --output_file ${EXP_ID}_feature_matrix.txt --prev_feature_matrix $MATRIX1 --next_feature_matrix $MATRIX2



if [ ! -f ${EXP_ID}_feature_matrix.txt ]; then
    echo "Building feature matrix $EXP_ID failed :("
fi

