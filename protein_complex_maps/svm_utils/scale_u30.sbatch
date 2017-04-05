#!/bin/bash
#SBATCH -J scalefull1           # job name
#SBATCH -o scalefull.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1             # total number of mpi tasks requested
#SBATCH -p normal    # queue (partition) -- normal, development, etc.
#SBATCH -t 4:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=cmcwhite@tacc.utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

/home1/03491/cmcwhite/libsvm/svm-scale -s atobsc3_u30full_libsvm1.scale_parameters atobsc3_u30full_libsvm1 > atobsc3_u30full_libsvm1.scale.txt



