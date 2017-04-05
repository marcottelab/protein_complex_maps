#!/bin/bash
#SBATCH -J predg0078125           # job name
#SBATCH -o predg0078125.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 130              # total number of mpi tasks requested
#SBATCH -p normal    # queue (partition) -- normal, development, etc.
#SBATCH -t 20:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

### Put the name of your commands file here.
export CONTROL_FILE=predict_commandsc32_g0078125.txt


### Nothing needs to change below here
module load launcher;
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=.
cd $WORKDIR;
$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE;


