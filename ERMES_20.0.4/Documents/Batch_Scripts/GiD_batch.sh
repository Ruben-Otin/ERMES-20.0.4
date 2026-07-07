###################################################################################################
# - Example of batch script to run ERMES on Linux bash using GiD offscreen.
# - See file "GiD_offscreen.bch" for an example of how to run GiD in batch mode.
###################################################################################################

# Set paths
PROBLEM_NAME=GiDProblemName
PROBLEM_PATH=$HOME/Projects/$PROBLEM_NAME.gid

GID_PATH=$HOME/GiDx64/gid-15.0.4
ERMES_EXE=$GID_PATH/problemtypes/ERMES_20.0.4/ERMES.gid/execs/unix/ERMES_20.0.4_c7.5
MACROS_PATH=$PROBLEM_PATH/GiDMacros
SOLVER_PATH=$HOME/Projects/PETSc

cd $PROBLEM_PATH

# Mesh with GiD offscreen
$GID_PATH/gid_offscreen -b $MACROS_PATH/Mesh.bch -n

# Delete previous files
rm -f $PROBLEM_NAME*dat Matrix_* Vector_* *.flavia.res *.info

# Calculate static currents after setting ERMES with GiD offscreen
$GID_PATH/gid_offscreen -b $MACROS_PATH/Calculate_J.bch -n
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info

# Move results files to a different folder
mv $PROBLEM_NAME.post.res Currents/$PROBLEM_NAME-J.post.res
mv $PROBLEM_NAME.post.msh Currents/$PROBLEM_NAME-J.post.msh 

# Delete previous files
rm -f $PROBLEM_NAME*dat Matrix_* Vector_* *.flavia.res *.info

# Calculate fields
$GID_PATH/gid_offscreen -b $MACROS_PATH/Calculate_E.bch -n
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info

# Solve with external solver
python3 $SOLVER_PATH/ERMES2PETSc.py