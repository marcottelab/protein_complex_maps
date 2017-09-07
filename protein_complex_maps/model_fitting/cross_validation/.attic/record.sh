
########
#Divide libsvm1 (training set) into segments
make_leaveout_train_div.sh
#Divides whole training file into 6 segments. 
#  Each possible set of 5 segments is used to train, predicting on the 6th
train_leaveout_divide.py
# This script is called by make_leaveout_train_div.sh

########
# Make training command file and run it on TACC
create_training_commands.sh
# A script to create the commands file for submitting to TACC
run_training.sbatch
# TACC submit file. Runs x_train_COMMANDS.sh made by creating_training_commands


########
# After training is done...
# Make testing command file and run it on TACC
create_testing_commands.sh
# A script to create the commands file for submitting to TACC
run_testing.sbatch
# TACC submit file. Runs x_leaveout_COMMANDS.sh made by creating_testing_commands

########
# Make pr command file and run it on TACC
create_prcurves_commands.sh
# Make a command file for making prcurves 
prcurve.py
# This script is called by make_pr_curves.s
run_prcurves.sbatch
# Tacc submission scripts. run x_pr_COMMANDS.sh 
