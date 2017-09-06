#!/bin/bash
#SBATCH -J scaletest1           # job name
#SBATCH -o scaletest.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -p normal    # queue (partition) -- normal, development, etc.
#SBATCH -t 1:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=cmcwhite@tacc.utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes


### Put the name of your commands file here.
export CONTROL_FILE=scale_test.sh

### Nothing needs to change below here
module load launcher;
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=.
cd $WORKDIR;
$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE;


