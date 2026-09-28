###############################################################################
# - This script computes the plasma permittivity tensor and outputs it to 
#   the binary file "Tensor_K_CPGE.bin" in the format expected by ERMES.
###############################################################################

import numpy as np

###############################################################################
# Plasma parameters
###############################################################################

f   = 28e9     ;       # Problem frequency
Wf  = 2*np.pi*f;       # Problem frequency (angular)

eo  = 8.8542e-12;      # Vacuum electric permittivity
Qe  = 1.6022e-19;      # Electron charge
Me  = 9.1095e-31;      # Electron mass

uKg = 1.660540199e-27; # Atomic mass (Kg)
MH1 = 1.0078250 * uKg; # Hidrogen mass
MH2 = 2.0141018 * uKg; # Deuterium mass
MH3 = 3.0160493 * uKg; # Tritium mass

###############################################################################
# Plasma species
###############################################################################
   
# Number of species in the plasma   
NumSpecies = 4;
    
# Charge, mass, and relative density
Qs = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) );  
Ms = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) ); 
Nr = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) ); 

Qs[ 0 ] =-Qe; Ms[ 0 ] = Me ; Nr[ 0 ] = 1.0; # Electron
Qs[ 1 ] = Qe; Ms[ 1 ] = MH1; Nr[ 1 ] = 0.0; # Hidrogen
Qs[ 2 ] = Qe; Ms[ 2 ] = MH2; Nr[ 2 ] = 1.0; # Deuterium
Qs[ 3 ] = Qe; Ms[ 3 ] = MH3; Nr[ 3 ] = 0.0; # Tritium

###############################################################################
# Damping
###############################################################################

eDamping = 1e9; 

Wd = np.zeros( NumSpecies, dtype=np.dtype( np.complex128 ) ); 

Wd[ 0 ] = Wf + 1j * eDamping;
Wd[ 1 ] = Wf;
Wd[ 2 ] = Wf;
Wd[ 3 ] = Wf;

###############################################################################
# Stix tensor components
###############################################################################

def Get_Stix_Tensor( ne, modB ):
    
    R = 1.0 + 1j * 0.0; 
    L = 1.0 + 1j * 0.0; 
    P = 1.0 + 1j * 0.0;
    
    for i in range( NumSpecies ):
        
        W2_ps   = Qs[i] * Qs[i] * ne * Nr[i] / ( eo * Ms[i] );
        Omega_s = Qs[i] * modB / Ms[i];
        
        R -= W2_ps / ( Wf * ( Wd[i] + Omega_s ) );
        L -= W2_ps / ( Wf * ( Wd[i] - Omega_s ) );
        P -= W2_ps / ( Wf * ( Wd[i]           ) );
    
    S = 0.5 * ( R + L );
    D = 0.5 * ( R - L );
    
    K = np.zeros( (3,3), dtype=np.dtype( np.complex128 ) );
    
    K[0][0] =  S  ; K[0][1] =+1j*D; K[0][2] = 0.0;
    K[1][0] =-1j*D; K[1][1] =  S  ; K[1][2] = 0.0;
    K[2][0] = 0.0 ; K[2][1] = 0.0 ; K[2][2] = P  ;
    
    return K;

###############################################################################
# Plasma files
###############################################################################

ne_FileName = r'./Plasma/ne_3D.dat'
Be_FileName = r'./Plasma/Be_3D.dat'

# Read electron density file
ne = [];
with open( ne_FileName, "r" ) as f:
    for line in f: 
        lsplit = line.split();
        
        ne.append( float( lsplit[1] ) );  
        
ne = np.array( ne );

# Read B field file
modB = [];
with open( Be_FileName, "r" ) as f:
    for line in f:
        lsplit = line.split();
        
        Bx = float( lsplit[1] );
        By = float( lsplit[2] );
        Bz = float( lsplit[3] );
        
        modB.append( np.sqrt( Bx**2 + By**2 + Bz**2 ) );
        
modB = np.array( modB );

###############################################################################
# Write Tensor_K_CPGE.bin file
###############################################################################

Tensor_K_FileName = r'./Tensor_K_CPGE.bin'

Kfile = open( Tensor_K_FileName, 'wb' );

for i in range( ne.size ):
    
    cK = Get_Stix_Tensor( ne[i], modB[i] );
    
    np.array( cK, dtype='<c16' ).tofile( Kfile );
       
Kfile.close();

###############################################################################



       