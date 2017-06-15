#!/bin/bash
#SBATCH -J runbuild     # job name
#SBATCH -o runbuild.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p largemem     # queue (partition) -- normal, development, etc.
#SBATCH -N 1
#SBATCH -c 16 		# number of CPUs/task
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

module load python/2.7.12


echo "Starting build join"
#Three is since I can't tell how the other job '2' is progressing
bash ../../features/build_feature_matrix.sh atobsc4_features.txt atobsc4


echo "Feature matrix built!"








