
#Doesn't seem to be a big different in log likelihood score  with different n's



python ../../topics/topics.py --infile plants_euNOG_concat.csv --sep="," > record.txt




#So, one set of IDs per line
#That's not so bad
#--cluster_predictions
#claire@marccomp01:/project/kdrew/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_prey_hein_prey_revisitTrain/clusterone_agglomod_clustering$ head blake_bioplex_prey_hein_prey_revisitTrain_corum_train_allComplexesCore_trainSplit_noTestOverlap_psweep_clusterone_agglomod.ii94.reduced.txt 
#2648 55689 10474 57325 8850 26009 6871
#8233 90324
#54943 515
#22908 9276 26958
#441502 3024
#8573 64130 55327 8777
#57476 113201 65983 9813



#--gold_standard 
#claire@marccomp01:/project/kdrew/data/protein_complex_maps/corum_revisit$ head allComplexesCore_geneid_merged06_testSplit1_humanOnly_under30.txt 
#7014 2957 2958
#898 1874 5933
#4088 4089 4087
#3937 5336
#2072 4436 2067
#3621 5111 2033
#1478 1479 672 580 1477
#54344 8813 8818
#3352 1901
#10010 7186 7187



#testing produced clusters using complex comparison.py

#Find out what these formats are, make topics.py match them
#Use training set, not test
#python complex_comparison.py --cluster_predictions ~/data/protein_complex_maps/v35_features/corum_split/v35_revisit/revisitTrain_fromHopper/blake_bioplex_prey_hein_prey_revisitTrain/trim_subunits/blake_bioplex_prey_hein_prey_revisitTrain_corum_train_allComplexesCore_trainSplit_noTestOverlap_psweep7.ii149.clusterone_agglomod.ii94.reduced.trimThreshold_reduced.txt --gold_standard ~/data/protein_complex_maps/corum_revisit/allComplexesCore_geneid_merged06_testSplit1_humanOnly_under30.txt --remove_non_gold_standard_proteins --normalize_by_combinations --pseudocount 0.00001





#
















#








