###############################################################################
# - This script defines functions and objects used in IRBC_Plane_Wave.py.
###############################################################################

import numpy as np
import re
from pathlib import Path

###############################################################################

class PlaneWave:  
    freq = 0.0;                    # Frequency

    modE = 0.0; phaE = 0.0;        # Module and phase of E

    ex = 0.0; ey = 0.0; ez = 0.0;  # Polarization vector
    kx = 0.0; ky = 0.0; kz = 0.0;  # Wave vector
    
    sgm_r = 0.0; sgm_i = 0.0;      # Conductivity 
    epr_r = 1.0; epr_i = 0.0;      # Relative permittivity
    mur_r = 1.0; mur_i = 0.0;      # Relative permeability
    
###############################################################################

class Node:   
    ID = 0;   # Node ID

    X = 0.0;  # Node cartesian coordinate X
    Y = 0.0;  # Node cartesian coordinate Y
    Z = 0.0;  # Node cartesian coordinate Z
    
###############################################################################

def Get_Element_Nodes_List( NodesIds, nXs, nYs, nZs ):
    """
    - Returns a list of element nodes objects.
    """
    eNodes = [];
    
    for i in range( len( NodesIds ) - 1 ):   
        Ni    = Node();  
        Ni.ID = NodesIds[ i ];
        
        Ni.X = nXs[ Ni.ID - 1 ];
        Ni.Y = nYs[ Ni.ID - 1 ];
        Ni.Z = nZs[ Ni.ID - 1 ];
        
        eNodes.append( Ni );
        
    return eNodes;

###############################################################################

def Read_Nodes_File():
    """
    - Reads nodes Ids and nodes coordinates from *-1.dat file.
    - Returns vectors: nId, nXc, nYc, nZc.
    """
    NodesFile = next( Path("..").glob("*-1.dat"), None );
    
    nId, nXc, nYc, nZc = [], [], [], [];
    
    pat = re.compile(r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;");
    
    with open( NodesFile, 'r' ) as f:
        for line in f:
            m = pat.search( line );
            if m:
                nId.append( int  ( m.group(1) ) );
                nXc.append( float( m.group(2) ) );
                nYc.append( float( m.group(3) ) );
                nZc.append( float( m.group(4) ) );
                
    nId = np.asarray( nId, dtype=np.int64 );
    nXc = np.asarray( nXc, float );
    nYc = np.asarray( nYc, float );
    nZc = np.asarray( nZc, float );

    return nId, nXc, nYc, nZc;

###############################################################################

def Read_IRBC_Elemets_File():
    """
    - Reads IRBC surface elements from *-19.dat file.
    - Returns element list: [ Node1, Node2, Node3, SurfaceId ].
    """
    ElementsFile = next( Path("..").glob("*-19.dat"), None );
    
    ElementList = [];
    
    with open( ElementsFile, 'r' ) as f:
        for line in f:
            iIds = list( map( int, line.split() ) );
            
            ElementList.append( [ iIds[0], iIds[1], iIds[2], iIds[3] ] );
  
    ElementList = np.asarray( ElementList, dtype=np.int64 );
    
    return ElementList;

###############################################################################

def Calculate_jK( PWave ):
    """
    - Returns jk = 1j*w*sqrt(ep*mu).
    """
    wf = 2.0 * np.pi * PWave.freq;
    w2 = wf * wf;

    eo = 8.8541878176e-12;
    mo = np.pi * 4.0e-7;            
    
    sg_real = PWave.sgm_r;
    sg_imag = PWave.sgm_i;

    ep_real = PWave.epr_r * eo; 
    ep_imag = PWave.epr_i * eo; 

    mu_real = PWave.mur_r * mo; 
    mu_imag = PWave.mur_i * mo; 

    cMu   = mu_real + 1j * mu_imag;
    cEpw2 = ep_real * w2 - sg_imag * wf + 1j * ( ep_imag * w2 + sg_real * wf );
    
    return ( 1j * np.sqrt( cEpw2 * cMu ) );    

###############################################################################

def Calculate_jK_div_Mu( PWave ):
    """
    - Returns jk/mu = 1j*w*sqrt(ep/mu).
    """
    wf = 2.0 * np.pi * PWave.freq;
    w2 = wf * wf;

    eo = 8.8541878176e-12;
    mo = np.pi * 4.0e-7;            
    
    sg_real = PWave.sgm_r;
    sg_imag = PWave.sgm_i;

    ep_real = PWave.epr_r * eo; 
    ep_imag = PWave.epr_i * eo; 

    mu_real = PWave.mur_r * mo; 
    mu_imag = PWave.mur_i * mo; 

    cMu   = mu_real + 1j * mu_imag;
    cEpw2 = ep_real * w2 - sg_imag * wf + 1j * ( ep_imag * w2 + sg_real * wf );
    
    return ( 1j * np.sqrt( cEpw2 / cMu ) );    

###############################################################################

def Calculate_Exterior_Normal( Nodes ):
    """
    - Returns the exterior normal of a surface element:  
      n = v2 x v1 / ||v2 x v1||
    """
    v1 = np.zeros( 3, dtype=np.dtype( np.float64 ) );
    v2 = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    v2[ 0 ] = Nodes[ 2 ].X - Nodes[ 0 ].X;
    v2[ 1 ] = Nodes[ 2 ].Y - Nodes[ 0 ].Y;
    v2[ 2 ] = Nodes[ 2 ].Z - Nodes[ 0 ].Z;
		    
    v1[ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X;
    v1[ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y;
    v1[ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z;

    Next = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
    Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
    Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

    Next = Next / np.linalg.norm( Next );
    
    return Next;

###############################################################################

def Calculate_Area( Nodes ):
    """
    - Returns the area of a surface element: 
      Area = 0.5 * ||v2 x v1||
    """
    v1 = np.zeros( 3, dtype=np.dtype( np.float64 ) );
    v2 = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    v2[ 0 ] = Nodes[ 2 ].X - Nodes[ 0 ].X;
    v2[ 1 ] = Nodes[ 2 ].Y - Nodes[ 0 ].Y;
    v2[ 2 ] = Nodes[ 2 ].Z - Nodes[ 0 ].Z;
		    
    v1[ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X;
    v1[ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y;
    v1[ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z;

    # Area = 0.5 * ||v2 x v1||
    Next = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
    Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
    Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

    Area = 0.5 * np.linalg.norm( Next );
    
    return Area;

###############################################################################

def Calculate_NiNj_Matrix( Nodes ):
    """
    - Returns Surf_Integral( NiNj ) matrix
    """
    NiNj_Matrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );
    
    Area = Calculate_Area( Nodes );
    
    for i in range( 3 ):
    	for j in range( 3 ):
            if( i != j ):
                NiNj_Matrix[ i ][ j ] = Area / 12.0;
            else:
                NiNj_Matrix[ i ][ j ] = Area / 6.0;		         

    return NiNj_Matrix;

###############################################################################

def Get_U_Vector_PW_RME_1st( PWave, Nodes ):
    """
    - Returns the residual vector for RME_1st element type.
    """
    # Nodes and Degrees of Freedom (DoF) of the element
    NumNodes = 3;
    NumDOFs  = NumNodes*3;
    
    # Initializing U vector
    U_vector = np.zeros( NumDOFs, dtype=np.dtype( np.complex128 ) );
    
    # Module and phase of E
    modul_F = PWave.modE; 
    phase_F = PWave.phaE;
    
    if( PWave.modE == 0.0 ): return U_vector;

    # E polarization vector ( normalized )
    Fpol = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Fpol[0] = PWave.ex; 
    Fpol[1] = PWave.ey; 
    Fpol[2] = PWave.ez;

    Fpol = Fpol / np.linalg.norm( Fpol );

    # Complex E vector
    cF = np.zeros( 3, dtype=np.dtype( np.complex128 ) );

    cF[0] = Fpol[0]; 
    cF[1] = Fpol[1]; 
    cF[2] = Fpol[2];

    cF = cF*( modul_F*np.cos( phase_F ) + 1j*modul_F*np.sin( phase_F ) );
    
    # Wave propagation vector ( normalized )
    Kdir = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Kdir[0] = PWave.kx;
    Kdir[1] = PWave.ky;
    Kdir[2] = PWave.kz;

    Kdir = Kdir / np.linalg.norm( Kdir );
    
    # Complex value jK
    jK = Calculate_jK( PWave );
    
    # Complex E vector on element nodes
    cFn = np.zeros( (NumNodes, NumNodes), dtype=np.dtype( np.complex128 ) );	
    
    for i in range( NumNodes ):
        KdirR  = Kdir[0]*Nodes[i].X + Kdir[1]*Nodes[i].Y + Kdir[2]*Nodes[i].Z;
        cFn[i] = cF;
        cFn[i] = cFn[i]*np.exp( jK * KdirR );
	                                                                         
	# Element exterior normal
    Next = Calculate_Exterior_Normal( Nodes );

	# k*n
    KdirNext = Kdir[0]*Next[0] + Kdir[1]*Next[1] + Kdir[2]*Next[2];

	# Complex value jK/Mu
    jK_div_Mu = Calculate_jK_div_Mu( PWave );
    
    # NiNj surface integrals
    NiNj_Matrix = Calculate_NiNj_Matrix( Nodes );
    
	# Curl element coefficient matrix ( n x n x E - n x k x E )
    CMatrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );

    CMatrix[0][0] = ( Next[0]*Next[0] - 1.0 ) - ( Kdir[0]*Next[0] - KdirNext );
    CMatrix[0][1] = ( Next[0]*Next[1]       ) - ( Kdir[0]*Next[1]            );
    CMatrix[0][2] = ( Next[0]*Next[2]       ) - ( Kdir[0]*Next[2]            );
    										  
    CMatrix[1][0] = ( Next[1]*Next[0]       ) - ( Kdir[1]*Next[0]            );
    CMatrix[1][1] = ( Next[1]*Next[1] - 1.0 ) - ( Kdir[1]*Next[1] - KdirNext );
    CMatrix[1][2] = ( Next[1]*Next[2]       ) - ( Kdir[1]*Next[2]            );
    										  
    CMatrix[2][0] = ( Next[2]*Next[0]       ) - ( Kdir[2]*Next[0]            );
    CMatrix[2][1] = ( Next[2]*Next[1]       ) - ( Kdir[2]*Next[1]            );
    CMatrix[2][2] = ( Next[2]*Next[2] - 1.0 ) - ( Kdir[2]*Next[2] - KdirNext );
    
    # Residual vector
    for i in range( NumNodes ):
        intg_NiCx = 0.0 + 1j * 0.0;
        intg_NiCy = 0.0 + 1j * 0.0;
        intg_NiCz = 0.0 + 1j * 0.0;

        for j in range( NumNodes ):
            intg_NiCx += ( CMatrix[0][0]*cFn[j][0] + 
                           CMatrix[0][1]*cFn[j][1] + 
                           CMatrix[0][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiCy += ( CMatrix[1][0]*cFn[j][0] + 
                           CMatrix[1][1]*cFn[j][1] + 
                           CMatrix[1][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiCz += ( CMatrix[2][0]*cFn[j][0] + 
                           CMatrix[2][1]*cFn[j][1] + 
                           CMatrix[2][2]*cFn[j][2] ) * NiNj_Matrix[i][j];

        U_vector[i             ] += jK_div_Mu * intg_NiCx;
        U_vector[i +   NumNodes] += jK_div_Mu * intg_NiCy;
        U_vector[i + 2*NumNodes] += jK_div_Mu * intg_NiCz;
        
    # Div element coefficient matrix ( n . E - k . E )
    DMatrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );

    DMatrix[0][0] = ( Next[0]*Next[0] ) - ( Next[0]*Kdir[0] );
    DMatrix[0][1] = ( Next[0]*Next[1] ) - ( Next[0]*Kdir[1] );
    DMatrix[0][2] = ( Next[0]*Next[2] ) - ( Next[0]*Kdir[2] );
								        
    DMatrix[1][0] = ( Next[1]*Next[0] ) - ( Next[1]*Kdir[0] );
    DMatrix[1][1] = ( Next[1]*Next[1] ) - ( Next[1]*Kdir[1] );
    DMatrix[1][2] = ( Next[1]*Next[2] ) - ( Next[1]*Kdir[2] );
									         
    DMatrix[2][0] = ( Next[2]*Next[0] ) - ( Next[2]*Kdir[0] );
    DMatrix[2][1] = ( Next[2]*Next[1] ) - ( Next[2]*Kdir[1] );
    DMatrix[2][2] = ( Next[2]*Next[2] ) - ( Next[2]*Kdir[2] );

    # Residual vector
    for i in range( NumNodes ):
        intg_NiDx = 0.0 + 1j * 0.0;
        intg_NiDy = 0.0 + 1j * 0.0;
        intg_NiDz = 0.0 + 1j * 0.0;

        for j in range( NumNodes ):
            intg_NiDx += ( DMatrix[0][0]*cFn[j][0] + 
                           DMatrix[0][1]*cFn[j][1] + 
                           DMatrix[0][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiDy += ( DMatrix[1][0]*cFn[j][0] + 
                           DMatrix[1][1]*cFn[j][1] + 
                           DMatrix[1][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiDz += ( DMatrix[2][0]*cFn[j][0] + 
                           DMatrix[2][1]*cFn[j][1] + 
                           DMatrix[2][2]*cFn[j][2] ) * NiNj_Matrix[i][j];

        U_vector[i             ] -= jK_div_Mu * intg_NiDx;
        U_vector[i +   NumNodes] -= jK_div_Mu * intg_NiDy;
        U_vector[i + 2*NumNodes] -= jK_div_Mu * intg_NiDz;

    return U_vector;

###############################################################################

def Get_P_Matrix_PW_RME_1st( PWave, Nodes ):
    """
    - Returns the stiffness matrix for RME_1st element type.
    """
    # Nodes and Degrees of Freedom (DoF) of the element
    NumNodes = 3;
    NumDOFs  = NumNodes*3;
    
	# Setting P matrix
    P_matrix = np.zeros((NumDOFs, NumDOFs), dtype=np.dtype(np.complex128));
    
    # Robin coefficient
    cRCoeff = -Calculate_jK_div_Mu( PWave );

    # NiNj surface integrals matrix
    NiNj_Matrix = Calculate_NiNj_Matrix( Nodes );

    # Stiffness matrix
    for i in range( NumNodes ):
        for j in range( NumNodes ):
            P_matrix[i             ][j             ] = cRCoeff*NiNj_Matrix[i][j];
            P_matrix[i +   NumNodes][j +   NumNodes] = cRCoeff*NiNj_Matrix[i][j];
            P_matrix[i + 2*NumNodes][j + 2*NumNodes] = cRCoeff*NiNj_Matrix[i][j];
            
    return P_matrix;

###############################################################################