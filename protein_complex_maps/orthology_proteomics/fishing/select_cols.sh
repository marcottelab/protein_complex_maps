#!/bin/bash
#SBATCH -J getcol      # job name
#SBATCH -o getcol.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p development     # queue (partition) -- normal, development, etc.
#SBATCH -c 16 		# number of CPUs/task
#SBATCH -t 2:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against


#Single columns already selected. This is for multiple columns...
awk -F',' '{print $1, $89, $90}' ../../features/atobsc_euNOG_corumtrain_labeled.txt > atobsc_euNOG_plants_poisson.txt




