************************************************************************************
* ERMES 20.0 - ERMES-PETSc direct interface
************************************************************************************

ERMESPETScSolver reads the ERMES linear system A Xo = B directly and in parallel,
solves it with PETSc, and writes the solution Vector_Xo.bin for ERMES. No Python and
no intermediate files are needed. Files read from the problem folder:

   Matrix_A_cmplx.bin, Matrix_A_int.bin          system matrix (coefficients, indices)
   Matrix_A_aux_cmplx.bin, Matrix_A_aux_int.bin  auxiliary matrix (Hermitic formats)
   Vector_B.bin                                  right-hand side

MATRIX FORMATS (detected automatically):

   Symmetric           Matrix_A: diagonal + upper diagonal of a complex symmetric matrix
   Full-matrix         Matrix_A: entire matrix
   Hermitic-Symmetric  Matrix_A: diagonal + upper diagonal of the Hermitian volumetric
                       contribution; Matrix_A_aux: diagonal + upper diagonal of the
                       complex symmetric Robin boundary contribution
   Hermitic-Full       Matrix_A: as above; Matrix_A_aux: entire Robin contribution,
                       added as stored

   Matrix_A_aux files present -> Hermitic format; the number of stored triangles selects
   the variant. The detected format is printed in the "*.info" file.

FILES:

   ERMESPETScSolver.cpp          Solver source code
   makefile                      Builds the solver
   ERMES2PETSc.sh                Script called by ERMES (solver settings at the top)
   win2wsl.sh                    Same, from Windows through WSL2

INSTALLATION:

1) Download PETSc from its git repository to the local folder "/petsc": 

   >> git clone -b release https://gitlab.com/petsc/petsc.git petsc

2) Configure PETSc using the command "./configure" inside the folder "/petsc" (check
   PETSc manual for configuration parameters definitions). Note that multiple 
   configurations can be installed on the same machine, for instance:

   >> ./configure PETSC_ARCH=arch-complex-M1 --with-debugging=0 --download-mumps 
      --download-scalapack --download-parmetis --download-metis --download-ptscotch 
      --with-64-bit-indices=1 --with-scalar-type=complex --download-mpich 
      --download-cmake --with-openmp --download-hwloc --with-cc=gcc --with-cxx=g++ 
      --with-fc=gfortran --download-fblaslapack --download-bison --download-make

3) Set PETSC_DIR and PETSC_ARCH before compiling and run PETSc (console or .bashrc):

   >> export PETSC_DIR=$HOME/petsc
   >> export PETSC_ARCH=arch-complex-M1

4) Compile the solver in this folder:

   >> make ERMESPETScSolver
   
5) Give execution permissions to "*.sh" files and "ERMESPETScSolver".

   >> chmod +x *.sh ERMESPETScSolver

USE FROM ERMES:

   Select "External solver" on "Solver type" and type on "External solvers settings":

   Linux  : bash /path/to/PETSc_Direct/ERMES2PETSc.sh
   Windows: wsl ../PETSc_Direct/win2wsl.sh

   Edit the settings at the top of ERMES2PETSc.sh (solver options). Inside a SLURM 
   job the solver uses all the nodes and tasks of the job. All the solver output, 
   including errors, is written to the ERMES "*.info" file.

USE WITHOUT ERMES:
   
   Define export and folder paths in a slurm.sh file: 
   
   >> export PETSC_DIR=$HOME/petsc
   >> export PETSC_ARCH=arch-complex-M1
   >> export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
   >> export OMP_NUM_THREADS=1
   >>
   >> SOLVER=$HOME/Scripts/PETSc_Direct/ERMESPETScSolver
   >> CASE_DIR=$HOME/projects/project.gid
   
   then use, for instance:

   >> $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 8 ./ERMESPETScSolver \
      -ermes_folder /path/to/Problem.gid \  
      -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps \
	  -memory_view -log_view 
	  
   >> $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n $SLURM_NTASKS $SOLVER \
      -ermes_folder $CASE_DIR \
      -ksp_type lgmres -pc_type sor -ksp_rtol 1e-6 -ksp_max_it 1000000 \
      -ermes_monitor_every 1000 -memory_view

   Solver options:
   -ermes_folder <dir>     : folder with the ERMES matrices and vector files.
   -ermes_monitor_every <n>: show residual every n iterations (iterative solvers).
   
   Useful MUMPS options:
   -pc_type cholesky     : save memory but only works for symmetric matrices.
   -mat_mumps_icntl_4 2  : print MUMPS info, including memory estimates.
   -mat_mumps_icntl_14 50: 50% more workspace (if MUMPS stops with INFOG(1)=-9).
   -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2: parallel ordering with ParMETIS.
   
NOTE:

   PETSc configured with --download-mpich must be launched with PETSc's own mpiexec
   (as in these scripts), not with srun.
   
   Check PETSc manual on "/External_Solvers/PETSc" for more info and solver options.
