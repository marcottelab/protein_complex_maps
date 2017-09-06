#!/bin/bash
#SBATCH -J a3ppformat           # job name
#SBATCH -o a3ppformat.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 100              # total number of mpi tasks requested
#SBATCH -p normal    # queue (partition) -- normal, development, etc.
#SBATCH -t 2:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

### This template is for running commands in parallel on TACC
### The commands file has one command per line

### Put the name of your commands file here.
export CONTROL_FILE=atobsc3_full_libsvm1.scale.poissonplus.formatting_COMMANDS.sh

### Nothing needs to change below here
module load launcher;
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=.
cd $WORKDIR;
$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE;
