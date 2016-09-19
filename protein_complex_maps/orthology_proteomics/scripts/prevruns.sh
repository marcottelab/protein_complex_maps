
#Plant map V2.0
#claire: running convert correlation matrix to pairs manually for broccoli nuclei
nohup python /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/convert_correlation.py --input_correlation_matrix Broccolinuclei_6-2016_raw_wide_elution_braol_euNOG.txt.corr_poisson --input_elution_profile Broccolinuclei_6-2016_raw_wide_elution_braol_euNOG.txt --output_file Broccolinuclei_6-2016_raw_wide_elution_braol_euNOG.txt.corr_poisson.pairs > nohup3.out&

#claire: running convert correlation matrix to pairs manually for selaginella
nohup python /project/cmcwhite/protein_complex_maps/protein_complex_maps/features/convert_correlation.py --input_correlation_matrix OP_SelaginellaSEC_20160309_raw_wide_elution_selml_euNOG.txt.corr_poisson --input_elution_profile OP_SelaginellaSEC_20160309_raw_wide_elution_selml_euNOG.txt --output_file OP_SelaginellaSEC_20160309_raw_wide_elution_selml_euNOG.txt.corr_poisson.pairs > nohup4.out&



