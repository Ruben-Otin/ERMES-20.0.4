###################################################################################################
# - Example of batch script to run ERMES on Linux bash.
###################################################################################################

# Set paths
ERMES_EXE=$HOME/GiDx64/gid-15.0.4/problemtypes/ERMES_20.0.4/ERMES.gid/execs/unix/ERMES_20.0.4_c7.5
PROBLEM_NAME=GiDProblemName
PROBLEM_PATH=$HOME/Projects/$PROBLEM_NAME.gid
SOLVER_PATH=$HOME/Projects/PETSc

# Set working directory
cd $PROBLEM_PATH

# Execute ERMES
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info

# Execute a python script calling PETSc
python3 $SOLVER_PATH/ERMES2PETSc.py