###############################################################################
# - This script calculates the stiffness matrix and residual vector of
#   every element on an imported Robin boundary condition surface for
#   the RME_1st element type.
#
# - The Robin boundary condition calculated on this script is a plane
#   wave boundary condition. 
#
# - The calculated stiffness matrices and residual vectors are saved 
#   in the binary files: "Matrix_P_IRBC.bin" and "Vector_U_IRBC.bin".
###############################################################################

import IRBC_Utils as iu
import numpy      as np

###############################################################################
# Plane wave properties.
###############################################################################

PWave = iu.PlaneWave();

# Frequency
PWave.freq  = 4e8;

# Module and phase of E
PWave.modE = 1.0; PWave.phaE = 0.0;

# Polarization vector
PWave.ex = 0.0; PWave.ey = 0.0; PWave.ez = 1.0;

# Wave vector
PWave.kx = 1.0; PWave.ky = 0.0; PWave.kz = 0.0;

# Material properties ( jk = j*w*sqrt(ep*mu) )
PWave.sgm_r = 0.0; PWave.sgm_i = 0.0;
PWave.epr_r = 1.0; PWave.epr_i = 0.0;
PWave.mur_r = 1.0; PWave.mur_i = 0.0;

###############################################################################
# Binary files paths and names.
###############################################################################

IRBC_U_FileName = r'../Vector_U_IRBC.bin'
IRBC_P_FileName = r'../Matrix_P_IRBC.bin'

###############################################################################
# Reading files, getting vector U and matrix P, and saving in files.
###############################################################################

# Read GiD mesh nodes:  No[ID] = p(x,y,z);
nIds, nXs, nYs, nZs = iu.Read_Nodes_File();

# Read IRBC surface elements
eList = iu.Read_IRBC_Elemets_File();

# Write in 'Vector_U_IRBC.bin' file
Ufile = open( IRBC_U_FileName, 'wb' );

# Write in 'Vector_P_IRBC.bin' file
Pfile = open( IRBC_P_FileName, 'wb' );

# Loop over IRBC surface elements
for i in range(len( eList )):
    
    # Get element nodes list
    eNodes = iu.Get_Element_Nodes_List( eList[i], nXs, nYs, nZs );
        
    # Get U vector and P matrix
    cU = iu.Get_U_Vector_PW_RME_1st( PWave, eNodes );
    cP = iu.Get_P_Matrix_PW_RME_1st( PWave, eNodes );
 
    # Write in binary files 
    np.array( cU, dtype='<c16' ).tofile( Ufile );
    np.array( cP, dtype='<c16' ).tofile( Pfile );
         
Ufile.close();
Pfile.close();

###############################################################################