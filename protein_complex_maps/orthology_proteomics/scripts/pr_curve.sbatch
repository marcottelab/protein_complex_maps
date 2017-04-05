#!/bin/bash
#SBATCH -J runsvm      # job name
#SBATCH -o runsvm.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -c 16 		# number of CPUs/task
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

module load python/2.7.12
# log2c=3 log2g=-1 rate=77.7885
BASEDIR=/scratch/03491/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics



python prcurve.py --results_wprob arathtraesorysj_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt arathtraesorysjbraolselml_euNOG_corum_train_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nodups_wprob.txt atobsc_euNOG_corumtrain_labeled.libsvm0.scaleByTrain.c8_g05_h0.resultsWprob_pairs_noself_nodups_wprob.txt --input_positives nonredundant_allComplexesCore_mammals_euNOG_merged06.test_ppis.txt --input_negatives nonredundant_allComplexesCore_mammals_euNOG_merged06.neg_test_ppis_12679.txt --output_file atobsc_euNOG_corumtrain_labeled.libsvm0.scaleByTrain.resultsWprob_pairs_noself_nondups_wprob.prcurve.pdf --labels arathtraesorysj arathtraesorysjbraolselml atobsc







