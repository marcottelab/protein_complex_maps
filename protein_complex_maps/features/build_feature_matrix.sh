
#nohup python -u /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files pvalcorrmin10_atobs_long_euNOG_concat_ln_pairs.txt  pvalcorrmin1_atobs_long_euNOG_concat_ln_pairs.txt  pvalcorrmin5_atobs_long_euNOG_concat_ln_pairs.txt  braol_euNOG_concat.txt.corr_poisson.pairs traes_euNOG_concat.txt.corr_poisson.pairs wgSEC1_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs  --output_file test_arathtraesorysjbraolselml_pval_euNOG_feature_matrix.txt >  nohup14.out&


#Put really large files last. 
#Speed up a lot!
#nohup python -u /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files pvalcorrmin10_atobs_long_euNOG_concat.txt_ln_pairs.txt  pvalcorrmin6_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin7_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin8_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin9_atobs_long_euNOG_concat_ln_pairs.txt At_Col_0_indark_201505_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_indark_fraction_201504_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_leaf_fraction_2014_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_leaf_fraction_2015_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs BroccoliNE_WWC_raw_wide_elution_braol_euNOG.txt.corr_poisson.pairs Broccolinuclei_6-2016_raw_wide_elution_braol_euNOG.txt.corr_poisson.pairs OP_SelaginellaSEC_20160309_raw_wide_elution_selml_euNOG.txt.corr_poisson.pairs orsyj_RiceL_IEX_raw_wide_elution_orysj_euNOG.csv_concat.txt.corr_poisson.pairs Rice_201505_uniprot_raw_wide_elution_orysj_euNOG.txt.corr_poisson.pairs RiceL_IEX_raw_wide_elution_orysj_euNOG.txt.corr_poisson.pairs selaginella_WWC_raw_wide_elution_selml_euNOG.txt.corr_poisson.pairs pvalcorrmin2_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin3_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin4_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin5_atobs_long_euNOG_concat_ln_pairs.txt arath_euNOG_concat.txt.corr_poisson.pairs braol_euNOG_concat.txt.corr_poisson.pairs selml_euNOG_concat.txt.corr_poisson.pairs wgIEF_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs wgSEC1_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs wheatgermIEX_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs WheatGermSEC_07-2015_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs traes_euNOG_concat.txt.corr_poisson.pairs pvalcorrmin1_atobs_long_euNOG_concat_ln_pairs.txt atobs_euNOG_concat.txt.corr_poisson.pairs --output_file atobs_pval_euNOG_feature_matrix.txt >  nohup16.out&
if [ $# -ne 2 ]; then
    echo usage: bash build_feature_matrix.sh new_experiment_basename features.txt
    echo ----new_experiment_basename: string that identifies this new feature matrix
    echo ----features.txt: file with 1 feature filename per line, sorted smallest file to largest
    echo ---------this makes building the feature matrix go faster
    exit 1
fi



FEATURES=$1
EXP_ID=$2
FEATURE_STR=$(cat $FEATURES | tr '\n' ' ')
echo ID: $EXP_ID
echo Features: $FEATURE_STR

nohup python -u  /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files $FEATURE_STR --output_file ${EXP_ID}_feature_matrix.txt >  nohup_${EXP_ID}_feature_matrix.out&



if [ ! -f ${EXP_ID}_feature_matrix.txt ]; then
    echo "Building feature matrix $EXP_ID with features $FEATURE_STR failed :("
fi


#What was run
#nohup python -u /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/build_feature_matrix.py --input_pairs_files pvalcorrmin10_atobs_long_euNOG_concat.txt_ln_pairs.txt pvalcorrmin1_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin2_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin3_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin4_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin5_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin6_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin7_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin8_atobs_long_euNOG_concat_ln_pairs.txt pvalcorrmin9_atobs_long_euNOG_concat_ln_pairs.txt arath_euNOG_concat.txt.corr_poisson.pairs arathtraesorysjbraolselml_euNOG_concat.txt.corr_poisson.pairs At_Col_0_indark_201505_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_indark_fraction_201504_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_leaf_fraction_2014_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs At_Col_0_leaf_fraction_2015_raw_wide_elution_arath_euNOG.txt.corr_poisson.pairs braol_euNOG_concat.txt.corr_poisson.pairs BroccoliNE_WWC_raw_wide_elution_braol_euNOG.txt.corr_poisson.pairs Broccolinuclei_6-2016_raw_wide_elution_braol_euNOG.txt.corr_poisson.pairs OP_SelaginellaSEC_20160309_raw_wide_elution_selml_euNOG.txt.corr_poisson.pairs orsyj_RiceL_IEX_raw_wide_elution_orysj_euNOG.csv_concat.txt.corr_poisson.pairs Rice_201505_uniprot_raw_wide_elution_orysj_euNOG.txt.corr_poisson.pairs RiceL_IEX_raw_wide_elution_orysj_euNOG.txt.corr_poisson.pairs selaginella_WWC_raw_wide_elution_selml_euNOG.txt.corr_poisson.pairs selml_euNOG_concat.txt.corr_poisson.pairs traes_euNOG_concat.txt.corr_poisson.pairs wgIEF_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs wgSEC1_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs wheatgermIEX_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs WheatGermSEC_07-2015_raw_wide_elution_traes_euNOG.txt.corr_poisson.pairs --output_file arathtraesorysjbraolselml_pval_euNOG_feature_matrix.txt >  nohup15.out&
