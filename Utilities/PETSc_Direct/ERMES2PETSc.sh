#!/bin/bash
###################################################################################################
#
#  ERMES2PETSc.sh  -  ERMES 20.0 / PETSc interface
#
#  Solves the ERMES linear system with ERMESPETScSolver, which reads the ERMES files directly in
#  parallel and writes the solution Vector_Xo.bin for ERMES.
#
#  Usage from ERMES:
#    1) Select "External solver" on ERMES "Solver type".
#    2) Type on ERMES "External solvers settings":
#       >> bash     /path_to/PETSc_Direct/ERMES2PETSc.sh (Linux)
#       >> wsl bash /path_to/PETSc_Direct/ERMES2PETSc.sh (Windows)
#
#  ERMES runs this script from the problem folder. The matrix format (Symmetric, Full-matrix,
#  Hermitic-Symmetric or Hermitic-Full) is detected automatically. Inside a SLURM job, the solver
#  runs on all the nodes and tasks of the job. All the solver output, including errors, goes to
#  the ERMES "*.info" file.
#
###################################################################################################

# -------------------------------------  User settings  -------------------------------------------

# PETSc configuration
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex-M1

# Solver type and options (see PETSc documentation). Examples:
# -Direct solver   : -ksp_type preonly -pc_type lu       -pc_factor_mat_solver_type mumps
# -Direct symmetric: -ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_type mumps
# -Iterative solver: -ksp_type lgmres  -pc_type sor -ksp_rtol 1e-6 -ksp_max_it 1000000
# -Options         : -ksp_monitor_true_residual (show residual for every iteration)
#                    -ermes_monitor_every 1000  (show residual every 1000 iterations)
#                    -memory_view -log_view     (consumed memory and solver info)  
SolverType="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
SolverOptions="-memory_view"

# Number of MPI processes (inside a SLURM job: all the tasks of the job)
NumParallTasks="${SLURM_NTASKS:-10}"

# Optional hybrid MPI + OpenMP: MPI processes per node and threads per process.
# Leave RanksPerNode empty to run NumParallTasks processes with 1 thread each.
# Example for 3 nodes x 32 cores: RanksPerNode="8", ThreadsPerRank="4" -> 24 processes x 4 threads.
RanksPerNode=""
ThreadsPerRank="1"

# MPI launcher: "petsc" (PETSc's own mpiexec, required with --download-mpich) | "srun" (only if
# PETSc was configured with the cluster MPI)
Launcher="petsc"

# Solver executable (default: in the same folder as this script)
SolverFullPath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/ERMESPETScSolver"

# -------------------------------------------------------------------------------------------------

# Send error messages to the standard output, so they also appear in the ERMES *.info file
exec 2>&1

export LD_LIBRARY_PATH="$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH"
FolderPath="$(pwd)"

# Processes per node and OpenMP threads
export OMP_NUM_THREADS="$ThreadsPerRank"
PerNodeFlag=""
if [ -n "$RanksPerNode" ]; then
    NumParallTasks=$(( RanksPerNode * ${SLURM_JOB_NUM_NODES:-1} ))
    if [ "$Launcher" = "srun" ]; then
        PerNodeFlag="--ntasks-per-node=$RanksPerNode --cpus-per-task=$ThreadsPerRank"
    else
        PerNodeFlag="-ppn $RanksPerNode"
    fi
fi

# MPI launcher
if [ "$Launcher" = "srun" ]; then
    MpiExec="srun"
else
    MpiExec="$PETSC_DIR/$PETSC_ARCH/bin/mpiexec"
    [ -x "$MpiExec" ] || MpiExec="$PETSC_DIR/lib/petsc/bin/petscmpiexec"
fi

# Solver settings
echo "------------------------------------------------------------"
echo "- PETSc  solver : $SolverType"
echo "- Solver options: $SolverOptions"
echo "- Processors    : $NumParallTasks x $ThreadsPerRank thread(s)"

# Remove any old solution, so that a failed run is never taken as valid by ERMES
rm -f "$FolderPath/Vector_Xo.bin"

# Solve. OMP_NUM_THREADS is also passed explicitly to every MPI process, because under SLURM
# PETSc's mpiexec does not always forward the environment (the processes would then use all the
# cores of the node as OpenMP threads, slowing down the factorization)
OmpFlag=""
[ "$Launcher" = "petsc" ] && OmpFlag="-genv OMP_NUM_THREADS $ThreadsPerRank"
"$MpiExec" -n "$NumParallTasks" $PerNodeFlag $OmpFlag "$SolverFullPath" \
    -ermes_folder "$FolderPath" $SolverType $SolverOptions
Status=$?

# A missing solution means the solver failed
if [ $Status -ne 0 ] || [ ! -s "$FolderPath/Vector_Xo.bin" ]; then
    echo "ERROR: PETSc solver failed (exit code $Status). No solution written."
    exit 1
fi
exit 0
