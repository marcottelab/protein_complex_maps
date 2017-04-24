#!/bin/bash
#SBATCH -J gridfull          # job name
#SBATCH -o gridfull.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 2              # total number of mpi tasks requested
#SBATCH -p largemem    # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

EXP_ID=atobsc3_full
PROGRAM_DIR=/home1/03491/cmcwhite

#$PROGRAM_DIR/libsvm/svm-scale -s ${EXP_ID}.libsvm1.scale_parameters ${EXP_ID}_libsvm1 > ${EXP_ID}_libsvm1.scale.txt
echo "starting" 
$PROGRAM_DIR/libsvm/tools/grid.py ${EXP_ID}_libsvm1.scale.txt


