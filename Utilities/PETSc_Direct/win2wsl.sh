#!/bin/bash
###################################################################################################
#
#  win2wsl.sh  -  Runs ERMES2PETSc.sh from Windows through WSL2.
#
#  Usage from ERMES: type on "External solvers settings":  wsl ../PETSc_Direct/win2wsl.sh
#  If an end-of-line error appears, run in WSL:  dos2unix win2wsl.sh ERMES2PETSc.sh
#
###################################################################################################

export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex-M1

bash "$(dirname "${BASH_SOURCE[0]}")/ERMES2PETSc.sh"
