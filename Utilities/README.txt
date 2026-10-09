***************************************************************************************************
* ERMES 20.0.4 - Utilities
***************************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - October 2026
***************************************************************************************************

This folder contains auxiliary scripts that complement ERMES 20.0.4. They help to automate
simulations on Linux and HPC clusters, solve the ERMES matrices with external solvers, and
generate files that can be imported into ERMES (plasma profiles, current sources, Robin
boundary conditions, element matrices). It includes the following folders:

- "Batch_Scripts"     : Examples of batch files to run GiD, ERMES and PETSc on Linux/SLURM.
- "Current_J_Scripts" : Generation of imported plasma current sources (J) from EQDSK files.
- "IRBC_Scripts"      : Generation of plane wave Imported Robin Boundary Conditions (IRBC).
- "IVEM_Scripts"      : Generation of imported element stiffness matrices (IVEM).
- "Mesh_Refiner"      : Parallel uniform refinement of the ERMES mesh (.dat files) and run.
- "NumPy_Solvers"     : Examples of Python NumPy/SciPy external solvers for ERMES.
- "PETSc_Direct"      : Parallel C++ PETSc solver that reads the ERMES system files directly.
- "Plasma_Scripts"    : EQDSK readers, cold plasma files generation and problem automation.
- "Scotty_GiD"        : Scotty beam tracing / vacuum Gaussian beams exported as GiD geometry.

All the scripts are examples intended to be copied and adapted to each specific problem.
Paths, file names, problem names and physical parameters are defined at the beginning of
each script and must be edited before use. Python scripts require Python 3 with NumPy and
SciPy (and Matplotlib for "eqdsk_reader.py" and the "Scotty_GiD" plots). The "Scotty_GiD"
scripts that run beam tracing also require Scotty (https://github.com/beam-tracing/Scotty).

***************************************************************************************************
1-) "Batch_Scripts"
***************************************************************************************************

Examples of shell scripts for running ERMES without the GiD graphical interface on Linux
workstations and HPC clusters.

- "ERMES_batch.sh"    : Runs the ERMES executable on an existing GiD problem and then calls
                        a Python script (ERMES2PETSc.py) to solve the system with PETSc.

- "GiD_batch.sh"      : Full workflow with GiD offscreen: meshes the problem, computes the
                        static currents, stores them in a "Currents" folder, computes the
                        fields, and solves with an external solver. It calls GiD batch
                        (.bch) macro files (Mesh.bch, Calculate_J.bch, Calculate_E.bch).

- "GiD_offscreen.bch" : Example of a GiD batch file. It creates a geometry, loads the ERMES
                        problemtype, defines Gaussian beam Robin coefficients and a cold
                        plasma material, sets the problem data (external PETSc solver),
                        assigns boundary conditions, meshes, and writes the calculation
                        file. Can also be run in GiD with "Import -> Batch file..." (Ctrl-b).

- "Linux_modules.sh"  : Example of Linux module configuration (gcc, cmake, python) for
                        ERMES. Adapt it to the modules available in your system.

- "PETSc_batch.sh"    : Example of a parallel PETSc solver call with petscmpiexec, including
                        solver options (e.g. MUMPS direct solver) and monitoring options.

- "SLURM_batch.sh"    : Example of SLURM job script for HPC clusters. It sets the PETSc
                        environment and runs ERMES in loops over several frequencies,
                        overwriting the ERMES input .dat files and collecting the results
                        into separate folders.

***************************************************************************************************
2-) "Current_J_Scripts"
***************************************************************************************************

Scripts that read an EQDSK plasma equilibrium file and compute the plasma current density
from the equilibrium profiles (p' and FF'). The result is written as an imported J source
file in the "Export_J_Sources" folder of the GiD problem. They require the EQDSK file and
the ERMES mesh files of the problem (nodes "*-1.dat" and volume elements "*-4.dat").

- "eqdsk2ermes_elevol_Jp.py" : J sources defined per tetrahedral volume element. The region
                               with current can be limited by coordinates or by material Id.

- "eqdsk2ermes_nodal_Jp.py"  : J sources defined per mesh node.

Note: by default these scripts read the EQDSK sample file from "./EQDSK_Samples/". Copy that
folder from "Plasma_Scripts" or edit the "EQDSK_File" path.

***************************************************************************************************
3-) "IRBC_Scripts"
***************************************************************************************************

- "IRBC_Plane_Wave_v2.py" : Generates a plane wave Imported Robin Boundary Condition for
                            EDG_1st and RME_1st elements. The user defines the frequency,
                            amplitude, phase, polarization, wave vector and material
                            properties of the plane wave. The script reads the mesh nodes
                            ("*-1.dat") and the IRBC surface elements ("*-19.dat") of the GiD
                            problem and writes the binary files "Vector_U_IRBC.bin" and
                            "Matrix_P_IRBC.bin" to be imported by ERMES.

***************************************************************************************************
4-) "IVEM_Scripts"
***************************************************************************************************

- "IVEM_K_Matrix.py" : Computes the element stiffness matrix of each tetrahedral element for
                       the RME_1st element type, given the material properties (frequency,
                       conductivity, permittivity, permeability). It reads the mesh files
                       ("*-1.dat" and "*-4.dat") and writes the binary file
                       "Matrix_K_IVEM.bin" to be imported by ERMES.

- "IVEM_Utils.py"    : Functions and objects used by "IVEM_K_Matrix.py" (mesh file reader,
                       tetrahedron shape function derivatives and volume, RME_1st element
                       matrix).

***************************************************************************************************
5-) "Mesh_Refiner"
***************************************************************************************************

Uniform refinement of an ERMES mesh directly on the .dat files written by GiD, without
re-meshing in GiD. It is useful to obtain very fine meshes (beyond what GiD can generate
comfortably) or to run convergence studies on the same geometry.

- "mesh_refiner.py" : Reads the nodes ("*-1.dat") and tetrahedra ("*-4.dat") of the problem,
                      builds a conforming refined mesh and rewrites every mesh-dependent
                      .dat file in the same format (nodal conditions "*-2.dat" and "*-3.dat",
                      volume and current elements "*-4.dat" and "*-5.dat", boundary
                      triangles "*-6.dat" to "*-8.dat", "*-10.dat" and "*-19.dat", prisms
                      "*-9.dat"). Original nodes keep their numbers, child elements inherit
                      the material, condition Id, etc. of their parent, and tetrahedra and
                      triangles keep their orientation. The other .dat files (materials,
                      plasma, frequency, solver...) are not modified. Two refinement modes:
                        'levels' : each edge divided in L parts (tetrahedron -> L^3,
                                   triangle -> L^2) without degrading the element quality.
                        'ps'     : Powell-Sabin / Worsey-Farin split (tetrahedron -> 24,
                                   triangle -> 6).
                      Runs in parallel on Linux, WSL and Windows. Settings are defined in the
                      USER SETTINGS block or on the command line, e.g.:
                        >> python mesh_refiner.py ./MyProblem.gid --level 2 --np 16
                        >> python mesh_refiner.py --help
                      IMPORTANT: by default the files are refined IN PLACE and the refinement
                      is NOT idempotent (running it twice refines the refined mesh). Let GiD
                      rewrite the .dat files first, or use "--out <dir>" to keep the input.

- "Refine2ERMES.sh" : Bash script (Linux or WSL) that runs the full pipeline inside the GiD
                      problem folder: removes previous results, converts the .dat files to
                      Unix line endings (dos2unix), refines the mesh with "mesh_refiner.py"
                      and runs ERMES on the refined mesh (log in the ERMES "*.info" file).
                      A stamp file (".refined_stamp") prevents refining the same mesh twice;
                      the refinement is re-enabled every time GiD rewrites the .dat files.
                      Problem name/folder, refiner and ERMES paths, mode, level, number of
                      processes and FORCE_REFINE (0 = auto, 1 = always, 2 = ERMES only) are
                      set at the top, or as arguments / environment variables, e.g.:
                        >> ./Refine2ERMES.sh MyProblem 3
                        >> LEVEL=4 NPROC=8 ./Refine2ERMES.sh
                        >> ./Refine2ERMES.sh -h

***************************************************************************************************
6-) "NumPy_Solvers"
***************************************************************************************************

Examples of Python scripts used as external solvers of ERMES. They read the matrix and
vector written by ERMES ("Matrix_A_cmplx.bin", "Matrix_A_int.bin", "Vector_B.bin"), solve
the system, and write the solution to "Vector_Xo.bin".

To use them, in the ERMES GiD interface select "Solving Parameters -> Solver -> Solver type
= External solver" and set: Solver executable = python, Solver parameters = path to the
script. If the matrix is stored in full format (not symmetric), comment the line that fills
the lower diagonal.

- "BiCG.py"    : Iterative BiConjugate Gradient solver with diagonal preconditioner.
- "SuperLU.py" : Direct sparse LU solver (SuperLU).

***************************************************************************************************
7-) "PETSc_Direct"
***************************************************************************************************

Direct interface between ERMES and PETSc. The C++ solver "ERMESPETScSolver" reads the ERMES
linear system in parallel (MPI-IO) directly from the ERMES binary files, solves it with any
PETSc solver (e.g. MUMPS direct solver or iterative KSP solvers), and writes the solution
"Vector_Xo.bin" read back by ERMES. No Python and no intermediate files are needed, and no
MPI process ever holds the whole matrix.

It reads "Matrix_A_cmplx.bin", "Matrix_A_int.bin" and "Vector_B.bin" and, for the Hermitic
formats, "Matrix_A_aux_cmplx.bin" and "Matrix_A_aux_int.bin". The matrix storage format
(Symmetric, Full-matrix, Hermitic-Symmetric or Hermitic-Full) is detected automatically and
printed in the ERMES "*.info" file. Requires PETSc >= 3.19 configured with complex scalars
and double precision.

- "ERMESPETScSolver.cpp" : Source code of the parallel PETSc solver. Accepts any PETSc option
                           plus "-ermes_folder <dir>" (folder with the ERMES files) and
                           "-ermes_monitor_every <n>" (print residual every n iterations).

- "makefile"             : Builds "ERMESPETScSolver" with the PETSc installation given by
                           PETSC_DIR and PETSC_ARCH (">> make ERMESPETScSolver").

- "ERMES2PETSc.sh"       : Script called by ERMES as external solver, on Linux
                           ("bash /path_to/PETSc_Direct/ERMES2PETSc.sh") and on Windows
                           through WSL ("wsl bash /path_to/PETSc_Direct/ERMES2PETSc.sh").
                           The PETSc configuration, solver type and options, number of MPI
                           processes, optional hybrid MPI+OpenMP and MPI launcher are set at
                           the top. Inside a SLURM job it uses all the nodes and tasks of the
                           job. All output and errors go to the ERMES "*.info" file. It can
                           also be run directly from a problem folder.

- "README.txt"           : Installation and usage instructions: PETSc download, Linux
                           modules, environment variables, example of PETSc configuration, 
                           compilation, use from ERMES, and solver options.

***************************************************************************************************
8-) "Plasma_Scripts"
***************************************************************************************************

Scripts for tokamak plasma problems (e.g. electron cyclotron waves in MAST-U).

- "eqdsk_reader.py"             : Reads an EQDSK equilibrium file and plots the poloidal flux,
                                  normalized flux, poloidal magnetic field, current density
                                  and electron density.

- "eqdsk2ermes_plasma_files.py" : Reads an EQDSK file and normalized density profiles and
                                  generates the electron density (ne) and magnetic field (Be)
                                  files required by the ERMES cold plasma material. Requires
                                  the list of nodes of the problem ("*-1.dat").

- "permittivity_tensor.py"      : Computes the cold plasma (Stix) permittivity tensor from the
                                  ne and Be files for a multi-species plasma and writes it to
                                  the binary file "Tensor_K_CPGE.bin" in ERMES format.

- "problem_generator.py"        : Automates a full simulation: computes the Gaussian beam
                                  launch geometry and O-mode polarization from the
                                  equilibrium, writes a GiD batch file, runs GiD offscreen to
                                  build and mesh the model, generates the plasma files, and
                                  runs ERMES. Windows and Linux paths are provided.

- "utilities.py"                : Functions used by "problem_generator.py" (vector rotations,
                                  beam and plasma geometry, GiD batch file writer, O-mode
                                  polarization, equilibrium file readers, plasma files).

- "EQDSK_Samples"               : Sample input data for the scripts above:
                                  - "mast-u-sample.eqdsk" : MAST-U EQDSK equilibrium file
                                  - "rho_ne.txt"          : Normalized electron density profile
                                  - "rho_te.txt"          : Normalized electron temperature profile
                                  - "rho_flux.txt"        : Normalized flux coordinate

***************************************************************************************************
9-) "Scotty_GiD"
***************************************************************************************************

Scripts that compute the path and shape of an electron cyclotron Gaussian beam and export it
as geometry that GiD can import (STL triangle mesh and/or exact NURBS surfaces in IGES). The
beam geometry can then be used to build and refine the ERMES model around the beam (e.g. the
MAST-U Doppler back-scattering launch). The beam envelope is the 1/e amplitude contour of the
beam. Parameters are set in the "USER PARAMETERS" block of each script, or on the command line
with "--name value" ("--show-params" and "--help" list them).

- "eqdsk_to_scotty.py" : Converts an EQDSK equilibrium file plus electron density and
                         temperature profiles into the input files of the Scotty beam tracer
                         ("topfile<sfx>", "ne<sfx>.dat", "Te<sfx>.dat" and, optionally,
                         "topfile<sfx>.json"). It also indicates the Scotty options to use
                         with these files (find_B_method, density_fit_method, etc.).

- "scotty_to_gid.py"   : Runs a Scotty beam tracing simulation (or reads an existing
                         "scotty_output*.h5" file) and writes the beam envelope as STL and/or
                         NURBS (IGES). The beam can be extended through vacuum beyond its first
                         and last points. It can also write an "FEM envelope": a closed
                         polyhedron with planar faces enclosing the beam at a given distance,
                         to be used as the ERMES computational volume. Text files with the
                         central ray, the E-field and wave vector k along the ray and at both
                         ends, and a summary of the beam parameters are written too.

- "beam_to_gid.py"     : Same outputs as "scotty_to_gid.py" for a Gaussian beam propagating
                         in vacuum, without Scotty and without plasma. The user defines the
                         frequency, start point, direction (vector or poloidal/toroidal
                         angles), beam profile (waist or launch width and curvature, also
                         elliptical/astigmatic beams), polarization and propagation length.

- "Plasma"             : Sample input data for "eqdsk_to_scotty.py":
                         - "mast-u-sample.eqdsk" : MAST-U EQDSK equilibrium file
                         - "flux_sample.dat"     : Normalized flux coordinate
                         - "ne_sample.dat"       : Electron density profile [m^-3]
                         - "te_sample.dat"       : Electron temperature profile [eV]

***************************************************************************************************

See "ERMES_20.0.4_Manual.pdf" for a description of the ERMES input/output files, imported
conditions, external solvers and batch execution used by these scripts.

ERMES 20.0.4 is licensed under the open-source 2-clause BSD license. Publications resulting
from the use of this software must cite the following article, which describes the program:

R. Otin, "ERMES 20.0: Open-source finite element tool for computational electromagnetics in
the frequency domain", Computer Physics Communications, Vol. 310, 109521, 2025.

For more information about ERMES 20.0.4, please contact: ruben.otin.bcn@gmail.com
