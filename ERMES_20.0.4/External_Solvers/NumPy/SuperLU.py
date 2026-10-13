##################################################################################################################
# 
# - Example of use of Python NumPy as an external solver of ERMES. 
# 
# - Before using this example, install Python and add "python.exe" to the operative system PATH. 
# 
# - In the ERMES GiD interface select: "Solving Parameters" -> "Solver" -> "Solver type" == "External solver".
# 
# - In the external solver settings textbox write:
#
#   1) Solver executable = python
#   2) Solver parameters = {python_script_path}\SuperLU.py (or "/" in Linux )
#
# - Remember to comment the matrix lower diagonal filling line in case of using a full matrix storage format.
#
##################################################################################################################

# Required imports
import numpy as np
import scipy.sparse.linalg 
import scipy.sparse
import sys
import os

# GiD problem path
Folder = os.getcwd();

if( sys.platform.startswith( 'win' ) ): 
    Folder = Folder + '\\';
else: 
    Folder = Folder + '/';

# Reading matrix  
B_vector = np.array( np.fromfile( Folder + r'Vector_B.bin'      , dtype=np.dtype( ( np.float64, 2 ) ) ) );
A_matrix = np.array( np.fromfile( Folder + r'Matrix_A_cmplx.bin', dtype=np.dtype( ( np.float64, 2 ) ) ) );
A_indexs = np.array( np.fromfile( Folder + r'Matrix_A_int.bin'  , dtype=np.dtype( ( np.int32  , 2 ) ) ) );

# Adapting formats
B_vector = B_vector[:,0] + 1j*B_vector[:,1];
A_matrix = A_matrix[:,0] + 1j*A_matrix[:,1];

A_matrix = scipy.sparse.csr_matrix( ( A_matrix, ( A_indexs[:,0]-1, A_indexs[:,1]-1 ) ) ); 

# Filling matrix lower diagonal for matrices stored in symmetric format.
# Comment this line in case of using full matrix storage format.
A_matrix = A_matrix + A_matrix.transpose() - scipy.sparse.diags( A_matrix.diagonal() ); 

# SuperLU solver
SuperLU  = scipy.sparse.linalg.splu( scipy.sparse.csr_matrix.tocsc( A_matrix ) );
X_vector = SuperLU.solve( B_vector );

# Showing solver info
print( '------------------------------------------------------------'                                       );
print( '- Python Super LU solver'                                                                           );
print( '- Residual: %r' % ( np.linalg.norm( A_matrix * X_vector - B_vector ) / np.linalg.norm( B_vector ) ) );
print( '------------------------------------------------------------'                                       );

# Writing solution vector on file
NewFile = open( Folder + r'Vector_Xo.bin', "wb" );
NewFile.write( X_vector );
NewFile.close();


