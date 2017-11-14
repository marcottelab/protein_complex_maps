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

BASEDIR=/scratch/03491/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics

bash template_run_svm.sh atobsc_euNOG $BASEDIR/correlation_elutions/features.txt $BASEDIR/correlation_elutions/atobsc_hyper_feature_matrix.txt $BASEDIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.train_ppis.txt $BASEDIR/corum/nonredundant_allComplexesCore_mammals_euNOG_merged06.neg_train_ppis_12679.txt

bash template_run_train.sh
