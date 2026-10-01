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

INSTALLATION:

1) Download PETSc from its git repository to the local folder "/petsc": 

   >> git clone -b release https://gitlab.com/petsc/petsc.git petsc
   
2) Set PETSC_DIR, PETSC_ARCH, and LD_LIBRARY_PATH before configure PETSc 
   (they must also be set before compiling and executing the solver):

   >> export PETSC_DIR=$HOME/petsc 
   >> export PETSC_ARCH=arch-complex-M1 
   >> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib
   
3) Configure PETSc using the command "./configure" inside the folder "/petsc" (check
   PETSc manual for configuration parameters definitions). Note that multiple 
   configurations can be installed on the same machine, for instance:

   >> ./configure PETSC_ARCH=arch-complex-M1 
      --with-cc=gcc 
      --with-cxx=g++ 
      --with-fc=gfortran 
      --with-debugging=0 
      --with-scalar-type=complex 
      --with-64-bit-indices=1 
      --with-openmp 
      --download-mpich 
      --download-hwloc 
      --download-openblas 
      --download-scalapack 
      --download-mumps 
      --download-metis 
      --download-parmetis 
      --download-ptscotch 
      --download-cmake 
      --download-bison 
      --download-make

4) Set PETSC_DIR, PETSC_ARCH, and LD_LIBRARY_PATH and compile solver on this folder:

   >> make ERMESPETScSolver
   
5) Give execution permissions to "*.sh" files and "ERMESPETScSolver".

   >> chmod +x *.sh ERMESPETScSolver

USE FROM ERMES:

   Select "External solver" on "Solver type" and type on "External solvers settings":

   Linux  : bash     /path_to/PETSc_Direct/ERMES2PETSc.sh
   Windows: wsl bash /path_to/PETSc_Direct/ERMES2PETSc.sh

   Edit the settings at the top of ERMES2PETSc.sh (solver options). Inside a SLURM 
   job the solver uses all the nodes and tasks of the job. All the solver output, 
   including errors, is written to the ERMES "*.info" file.
   
   Solver options:
   -ermes_folder <dir>     : folder with the ERMES matrices and vector files.
   -ermes_monitor_every <n>: show residual every n iterations (iterative solvers).
   
   Useful MUMPS options:
   -pc_type cholesky     : save memory but only works for symmetric matrices.
   -mat_mumps_icntl_4 2  : print MUMPS info, including memory estimates.
   -mat_mumps_icntl_14 50: 50% more workspace (if MUMPS stops with INFOG(1)=-9).
   -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2: parallel ordering with ParMETIS.
   
   PETSc can also be called directly by executing ERMES2PETSc.sh and editing:
   - Solver parameters (solver type, number of tasks,etc.)
   - SolverFullPath    (path to "ERMESPETScSolver")
   - FolderPath        (path to the linear system files)
  
NOTE:

   PETSc configured with --download-mpich must be launched with PETSc's own mpiexec
   (as in these scripts), not with srun.
   
   Check PETSc manual on "/External_Solvers/PETSc" or visit https://petsc.org for 
   more info and solver options.
