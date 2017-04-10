#!/bin/bash
#SBATCH -J runsvm      # job name
#SBATCH -o runsvm.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p largemem     # queue (partition) -- normal, development, etc.
#SBATCH -c 16 		# number of CPUs/task
#SBATCH -t 3:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

#sed -E 's/"|\[|\(|\]|\)//g' atobsc_euNOG_corumtrain_labeled.tmp > atobsc_euNOG_corumtrain_labeled.tmp2

#sed -E "s/'//g" atobsc_euNOG_corumtrain_labeled.tmp2 > atobsc_euNOG_corumtrain_labeled.tmp3
#sed -E "s/, / /g" atobsc_euNOG_corumtrain_labeled.tmp3 > atobsc_euNOG_corumtrain_labeled.tmp4


module load python/2.7.12

#python ../../external/score.py plants_euNOG_concat.txt spearman

python ../../features/convert_correlation.py --input_correlation_matrix plants_euNOG_concat.txt.corr_spearman --input_elution_profile plants_euNOG_concat.txt --output_file  plants_euNOG_concat.txt.corr_spearman.pairs 


python ../../features/alphabetize_pairs.py --feature_pairs plants_euNOG_concat.txt.corr_spearman.pairs --outfile plants_euNOG_concat.txt.corr_spearman.pairs.ordered 


