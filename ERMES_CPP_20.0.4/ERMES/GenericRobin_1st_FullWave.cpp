
#include "GenericRobin_1st_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //**********************************************************************************************************************
    //*******************************************  Quasi-Static Mode  ******************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Real part of the generic Robin coefficient
    // COMPLEX_IBC[ 1 ] = Imag part of the generic Robin coefficient
    // 
    // COMPLEX_IBC_2o[ 0 ] = Flux module
    // COMPLEX_IBC_2o[ 1 ] = Flux phase
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 0
    //  
    //**********************************************************************************************************************
    //*******************************************  Plane Wave Mode  ********************************************************
    //**********************************************************************************************************************    
    //
    // COMPLEX_IBC[ 0 ] = Real part of the electric conductivity in S/m
    // COMPLEX_IBC[ 1 ] = Imag part of the electric conductivity in S/m
    // 
    // COMPLEX_IBC_2o[ 0 ] = Real part of the relative electric permittivity
    // COMPLEX_IBC_2o[ 1 ] = Imag part of the relative electric permittivity
    //  
    // COMPLEX_NEUMANN_FLOW[ 0 ] = Real part of the relative magnetic permeability
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Imag part of the relative magnetic permeability
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 1 in single   plane wave  mode  
    //                           = 2 in multiple plane waves mode 
    //                            
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector
    // 
    //**********************************************************************************************************************
    //*******************************************  Gaussian Beam Mode  *****************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Module of the E field perpendicular to the polarization and wave vectors at the waist center
    // COMPLEX_IBC[ 1 ] = Phase  of the E field perpendicular to the polarization and wave vectors at the waist center
    // 
    // COMPLEX_IBC_2o[ 0 ] = Gaussian beam waist radius on the direction parallel to the polarization vector
    // COMPLEX_IBC_2o[ 1 ] = Gaussian beam waist radius on the direction perpendicular to the polarization and wave vectors
    //   
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of the Gaussian beam waist 
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector at the waist center
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector at the waist center
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 3
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector 
    // 
    //**********************************************************************************************************************
    
    //**********************************************************************************************************************
    //* - Global IDs of the DOFs
    //**********************************************************************************************************************	
    void GenericRobin_1st_FullWave::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

		if( mPotentialsOn )
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcAz()->EquationId();
                EquationId[ i + 3*mNumNodes ] = mNodes[ i ]->pDofcVs()->EquationId();
            }				
		}
		else
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcEy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcEz()->EquationId();
            }		
		}
    }

    //**********************************************************************************************************************
    // - Area of the element
    //**********************************************************************************************************************
    double GenericRobin_1st_FullWave::Calculate_Area()
    {
		double na[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        // Area = 0.5 * ||v2 x v1||
        na[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        na[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        na[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

		return ( 0.5 * sqrt( na[ 0 ]*na[ 0 ] + na[ 1 ]*na[ 1 ] + na[ 2 ]*na[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - Exterior normal 
    //**********************************************************************************************************************	
    void GenericRobin_1st_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
    {
		double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        Next.resize( 3 );

        Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        Next *= ( 1.0 / sqrt( Next[ 0 ]*Next[ 0 ] + Next[ 1 ]*Next[ 1 ] + Next[ 2 ]*Next[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - jk / mu = j * w * sqrt ( ep / mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_1st_FullWave::Calculate_jK_div_Mu()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo / mo ) );
        }
        // Plane wave 
        else
        {
            if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] == 0.0 ) 
            {
                return std::complex<double>( 0.0, 0.0 );
            }            
            
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );

            return ( cUnit * sqrt( cEpw2 / cMu ) );
        }
    }

    //**********************************************************************************************************************
    //* - jk = j * w * sqrt ( ep * mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_1st_FullWave::Calculate_jK()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo * mo ) );
        }
        // Plane wave 
        else
        {
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );
            
            return ( cUnit * sqrt( cEpw2 * cMu ) );
        }
    }

    //**********************************************************************************************************************
    //* - NiNj surface integrals matrix
    //**********************************************************************************************************************	
    void GenericRobin_1st_FullWave::Calculate_NiNj_Matrix( Matrix<double>& NiNj_Matrix )
    {
        NiNj_Matrix.Resize( mNumNodes, mNumNodes );

        for( int i=0; i<mNumNodes; i++ )
        {
        	for( int j=0; j<mNumNodes; j++ )
        	{
        	    if( i != j )
                {
					NiNj_Matrix[ i ][ j ] = mArea / 12.0;
                }
        	    else
                {
					NiNj_Matrix[ i ][ j ] = mArea / 6.0;		         
                }
        	}
        }
	}

    //**********************************************************************************************************************
    //* - Residual vector and quasi-static flux values ( needed for total J sources calculations )
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetResidualVector_Flux( std::map< unsigned int, std::complex<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

        if( ( mPotentialsOn == true ) && ( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 ) ) 
		{
            double Intg_Ni = mArea / 3.0;
            
            Vector<double> Coef = (*mProperties)( COMPLEX_IBC    );
            Vector<double> Flux = (*mProperties)( COMPLEX_IBC_2o );

            std::complex<double> cFlux( Flux[ 0 ] * cos( Flux[ 1 ] ), Flux[ 0 ] * sin( Flux[ 1 ] ) );

            for( int i=0; i<3; i++ ) 
            { 
                if( ( Coef[ 0 ] == 0.0 ) && ( Coef[ 1 ] == 0.0 ) ) 
                {
                    FluxValue[ mNodes[ i ]->Id() ] = cFlux;
                }
                
                ResidualVector[ i + 3*mNumNodes ] += cFlux * Intg_Ni;
            }
		}
    }

    //**********************************************************************************************************************
    //* - Residual vector for a single plane wave
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetResidualVector_Single_Wave( Vector< std::complex<double> >& ResidualVector )
    {
	    // Module and phase of the primary field ( E/A )
        double modul_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double phase_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];

		if( modul_F == 0.0 ) return;

		// E/A polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Complex E/A vector
        Vector< std::complex<double> > cF( 3 );

        cF[ 0 ] = Fpol[ 0 ]; 
        cF[ 1 ] = Fpol[ 1 ]; 
        cF[ 2 ] = Fpol[ 2 ];

        cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) );

		if( mPotentialsOn ) 
		{
			cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

		// Complex value jK
        std::complex<double> jK = Calculate_jK();

		// Complex E/A vector on element nodes
		Vector< Vector< std::complex<double> > > cFn( mNumNodes );

		for( int i=0; i<mNumNodes; i++ )
		{
			double KdirR = Kdir[ 0 ]*mNodes[ i ]->X() + Kdir[ 1 ]*mNodes[ i ]->Y() + Kdir[ 2 ]*mNodes[ i ]->Z();
			
			cFn[ i ]  = cF;

			cFn[ i ] *= std::exp( jK * KdirR );
		}

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

        // NiNj surface integrals
        Matrix<double> NiNj_Matrix;

		Calculate_NiNj_Matrix( NiNj_Matrix );

		// Curl element coefficient matrix ( n x n x E - n x k x E )
        Matrix<double> CMatrix( 3, 3 );

		CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
										  
		CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
										  
		CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg_NiCx( 0.0, 0.0 );
			std::complex<double> intg_NiCy( 0.0, 0.0 );
			std::complex<double> intg_NiCz( 0.0, 0.0 );

            for( int j=0; j<mNumNodes; j++ )
			{
				intg_NiCx += ( CMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiCy += ( CMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiCz += ( CMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
			}

            ResidualVector[ i               ] += jK_div_Mu * intg_NiCx;
            ResidualVector[ i +   mNumNodes ] += jK_div_Mu * intg_NiCy;
			ResidualVector[ i + 2*mNumNodes ] += jK_div_Mu * intg_NiCz;
        }  

		// Div element coefficient matrix ( n . E - k . E )
        Matrix<double> DMatrix( 3, 3 );

        DMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] ) - ( Next[ 0 ]*Kdir[ 0 ] );
		DMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ] ) - ( Next[ 0 ]*Kdir[ 1 ] );
		DMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ] ) - ( Next[ 0 ]*Kdir[ 2 ] );
								        
		DMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ] ) - ( Next[ 1 ]*Kdir[ 0 ] );
		DMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] ) - ( Next[ 1 ]*Kdir[ 1 ] );
		DMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ] ) - ( Next[ 1 ]*Kdir[ 2 ] );
									         
		DMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ] ) - ( Next[ 2 ]*Kdir[ 0 ] );
		DMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ] ) - ( Next[ 2 ]*Kdir[ 1 ] );
		DMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] ) - ( Next[ 2 ]*Kdir[ 2 ] );

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg_NiDx( 0.0, 0.0 );
			std::complex<double> intg_NiDy( 0.0, 0.0 );
			std::complex<double> intg_NiDz( 0.0, 0.0 );

            for( int j=0; j<mNumNodes; j++ )
			{
				intg_NiDx += ( DMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiDy += ( DMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiDz += ( DMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
			}

            ResidualVector[ i               ] -= jK_div_Mu * intg_NiDx;
            ResidualVector[ i +   mNumNodes ] -= jK_div_Mu * intg_NiDy;
			ResidualVector[ i + 2*mNumNodes ] -= jK_div_Mu * intg_NiDz;
        }  
    }

    //**********************************************************************************************************************
    //* - Residual vector for multiple plane waves
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetResidualVector_Multpl_Waves( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
	{
        // If plane waves set is empty then return
		if( FluxValue.size() == 0 ) return; 

	    // Complex value jK
        std::complex<double> jK = Calculate_jK();

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

        // NiNj surface integrals matrix
        Matrix<double> NiNj_Matrix;
		
		Calculate_NiNj_Matrix( NiNj_Matrix );

		// Iterator over plane waves
		Vector< Vector<double> >::iterator pw_it;

        for( pw_it = FluxValue.begin(); pw_it != FluxValue.end(); ++pw_it )
		{
            // Module and phase of the primary field ( E/A )
            double modul_F = (*pw_it)[ 0 ]; 
            double phase_F = (*pw_it)[ 1 ];

			if( modul_F == 0.0 ) continue;
		    
		    // E/A polarization vector ( normalized )
            Vector<double> Fpol( 3 );
		    
            Fpol[ 0 ] = (*pw_it)[ 2 ];
            Fpol[ 1 ] = (*pw_it)[ 3 ];
            Fpol[ 2 ] = (*pw_it)[ 4 ];
		    
            Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );
		    
            // Complex E/A vector
            Vector< std::complex<double> > cF( 3 );
		    
            cF[ 0 ] = Fpol[ 0 ]; 
            cF[ 1 ] = Fpol[ 1 ]; 
            cF[ 2 ] = Fpol[ 2 ];
		    
            cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) );
		    
		    if( mPotentialsOn ) 
		    {
		    	cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		    }
		    
            // Wave propagation vector ( normalized )
            Vector<double> Kdir( 3 );
		    
            Kdir[ 0 ] = (*pw_it)[ 5 ];
            Kdir[ 1 ] = (*pw_it)[ 6 ];
            Kdir[ 2 ] = (*pw_it)[ 7 ];
		    
            Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );
		    
		    // Complex E/A vector on element nodes
		    Vector< Vector< std::complex<double> > > cFn( mNumNodes );
		    
		    for( int i=0; i<mNumNodes; i++ )
		    {
		    	double KdirR = Kdir[ 0 ]*mNodes[ i ]->X() + Kdir[ 1 ]*mNodes[ i ]->Y() + Kdir[ 2 ]*mNodes[ i ]->Z();
		    	
		    	cFn[ i ]  = cF;
		    
		    	cFn[ i ] *= std::exp( jK * KdirR );
		    }
		    
		    // k*n
		    double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];
		    
		    // Curl element coefficient matrix ( n x n x E - n x k x E )
            Matrix<double> CMatrix( 3, 3 );
		    
		    CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		    CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		    CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
		    		 	 					  
		    CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		    CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		    CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
		    		 	 					  
		    CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		    CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		    CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );
		    
	        // Residual vector
		    for( int i=0; i<mNumNodes; i++ )
            {
                std::complex<double> intg_NiCx( 0.0, 0.0 );
		    	std::complex<double> intg_NiCy( 0.0, 0.0 );
		    	std::complex<double> intg_NiCz( 0.0, 0.0 );
		    
                for( int j=0; j<mNumNodes; j++ )
		    	{
		    		intg_NiCx += ( CMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    		intg_NiCy += ( CMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    		intg_NiCz += ( CMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    	}
		    
                ResidualVector[ i               ] += jK_div_Mu * intg_NiCx;
                ResidualVector[ i +   mNumNodes ] += jK_div_Mu * intg_NiCy;
		    	ResidualVector[ i + 2*mNumNodes ] += jK_div_Mu * intg_NiCz;
            }  

		    // Div element coefficient matrix ( n . E - k . E )
            Matrix<double> DMatrix( 3, 3 );
            
            DMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] ) - ( Next[ 0 ]*Kdir[ 0 ] );
		    DMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ] ) - ( Next[ 0 ]*Kdir[ 1 ] );
		    DMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ] ) - ( Next[ 0 ]*Kdir[ 2 ] );
		    						        
		    DMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ] ) - ( Next[ 1 ]*Kdir[ 0 ] );
		    DMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] ) - ( Next[ 1 ]*Kdir[ 1 ] );
		    DMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ] ) - ( Next[ 1 ]*Kdir[ 2 ] );
		    							         
		    DMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ] ) - ( Next[ 2 ]*Kdir[ 0 ] );
		    DMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ] ) - ( Next[ 2 ]*Kdir[ 1 ] );
		    DMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] ) - ( Next[ 2 ]*Kdir[ 2 ] );
            
	        // Residual vector
		    for( int i=0; i<mNumNodes; i++ )
            {
                std::complex<double> intg_NiDx( 0.0, 0.0 );
		    	std::complex<double> intg_NiDy( 0.0, 0.0 );
		    	std::complex<double> intg_NiDz( 0.0, 0.0 );
            
                for( int j=0; j<mNumNodes; j++ )
		    	{
		    		intg_NiDx += ( DMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    		intg_NiDy += ( DMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    		intg_NiDz += ( DMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
		    	}
            
                ResidualVector[ i               ] -= jK_div_Mu * intg_NiDx;
                ResidualVector[ i +   mNumNodes ] -= jK_div_Mu * intg_NiDy;
		    	ResidualVector[ i + 2*mNumNodes ] -= jK_div_Mu * intg_NiDz;
            }  
        }
	}

    //**********************************************************************************************************************
    //* - Residual vector for Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetResidualVector_Gauss_Beams( Vector< std::complex<double> >& ResidualVector )
    {
        // Gaussian beam waist radius 
        double GBR_pol = (*mProperties)( COMPLEX_IBC_2o )[ 0 ]; 
        double GBR_per = (*mProperties)( COMPLEX_IBC_2o )[ 1 ];

		if( ( GBR_pol == 0.0 ) || ( GBR_per == 0.0 ) ) return;

		// Field polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

        // Fper = Kdir X Fpol ( normalized )
        Vector<double> Fper( 3 );

        Fper[ 0 ] = Kdir[ 1 ]*Fpol[ 2 ] - Kdir[ 2 ]*Fpol[ 1 ];
        Fper[ 1 ] = Kdir[ 2 ]*Fpol[ 0 ] - Kdir[ 0 ]*Fpol[ 2 ];
        Fper[ 2 ] = Kdir[ 0 ]*Fpol[ 1 ] - Kdir[ 1 ]*Fpol[ 0 ];

        Fper *= ( 1.0 / sqrt( Fper[ 0 ]*Fper[ 0 ] + Fper[ 1 ]*Fper[ 1 ] + Fper[ 2 ]*Fper[ 2 ] ) );

        // Field parallel to the polarization vector
        Vector< std::complex<double> > Fpar( 3 );

        Fpar[ 0 ] = Fpol[ 0 ]; 
        Fpar[ 1 ] = Fpol[ 1 ]; 
        Fpar[ 2 ] = Fpol[ 2 ];

        double mod_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double pha_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
        
        Fpar *= std::complex<double>( mod_Fpar * cos( pha_Fpar ), mod_Fpar * sin( pha_Fpar ) );

        // Field perpendicular to the polarization and wave vectors
        Vector< std::complex<double> > Fppw( 3 );

        Fppw[ 0 ] = Fper[ 0 ]; 
        Fppw[ 1 ] = Fper[ 1 ]; 
        Fppw[ 2 ] = Fper[ 2 ];

        double mod_Fppw = (*mProperties)( COMPLEX_IBC )[ 0 ]; 
        double pha_Fppw = (*mProperties)( COMPLEX_IBC )[ 1 ];
        
        Fppw *= std::complex<double>( mod_Fppw * cos( pha_Fppw ), mod_Fppw * sin( pha_Fppw ) );

        // Total field on waist central point
        Vector< std::complex<double> > Ftot( 3 );

        Ftot[ 0 ] = Fpar[ 0 ] + Fppw[ 0 ]; 
        Ftot[ 1 ] = Fpar[ 1 ] + Fppw[ 1 ]; 
        Ftot[ 2 ] = Fpar[ 2 ] + Fppw[ 2 ];

		if( mPotentialsOn ) 
		{
			Ftot *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

		// Complex wave vector jK 
        std::complex<double> jK = Calculate_jK();

        // Rayleight range
        double ZR_pol = 0.5 * std::imag( jK ) * GBR_pol * GBR_pol;
        double ZR_per = 0.5 * std::imag( jK ) * GBR_per * GBR_per;

        // Gaussian beam waist center coordinates
        double Xo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];			
        double Yo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
        double Zo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ]; 

		// Field vector on element nodes
		Vector< Vector< std::complex<double> > > cFn( mNumNodes );

		for( int i=0; i<mNumNodes; i++ )
		{
            // Node coordinates
			double Xn = mNodes[ i ]->X();
            double Yn = mNodes[ i ]->Y();
            double Zn = mNodes[ i ]->Z();

            // Node position from waist center
            Vector<double> Rn( 3 );

            Rn[ 0 ] = Xn - Xo; 
            Rn[ 1 ] = Yn - Yo; 
            Rn[ 2 ] = Zn - Zo;

            // Node distance from waist plane
            double Zd = Kdir[ 0 ]*Rn[ 0 ] + Kdir[ 1 ]*Rn[ 1 ] + Kdir[ 2 ]*Rn[ 2 ];

            // Accumulated phase 
            std::complex<double> AccPhase = jK * Zd;

            // Gouy phase shift
            std::complex<double> GouyPhase( 0.0, 0.5 * ( std::atan( Zd / ZR_pol ) + std::atan( Zd / ZR_per ) ) );

            // Q parameters
            std::complex<double> Qpol( Zd, -ZR_pol );
            std::complex<double> Qper( Zd, -ZR_per );

            // Node coordinates projection on waist plane
            double Xpol = Fpol[ 0 ]*Rn[ 0 ] + Fpol[ 1 ]*Rn[ 1 ] + Fpol[ 2 ]*Rn[ 2 ];
            double Xper = Fper[ 0 ]*Rn[ 0 ] + Fper[ 1 ]*Rn[ 1 ] + Fper[ 2 ]*Rn[ 2 ];

            // Q parameters phase
            std::complex<double> QpPhase = jK * 0.5 * ( ( Xpol * Xpol / Qpol ) + ( Xper * Xper / Qper ) );

            // Field normalization 
            double Beta = std::sqrt( std::sqrt( ZR_pol * ZR_per / ( std::abs( Qpol ) * std::abs( Qper ) ) ) ); 

            // Field value on node 
            cFn[ i ]  = Ftot;
            cFn[ i ] *= Beta * std::exp( AccPhase - GouyPhase + QpPhase );
		}

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

        // NiNj surface integrals
        Matrix<double> NiNj_Matrix;

		Calculate_NiNj_Matrix( NiNj_Matrix );

		// Curl element coefficient matrix ( n x n x E - n x k x E )
        Matrix<double> CMatrix( 3, 3 );

		CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
										  
		CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
										  
		CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg_NiCx( 0.0, 0.0 );
			std::complex<double> intg_NiCy( 0.0, 0.0 );
			std::complex<double> intg_NiCz( 0.0, 0.0 );

            for( int j=0; j<mNumNodes; j++ )
			{
				intg_NiCx += ( CMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiCy += ( CMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiCz += ( CMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + CMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + CMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
			}

            ResidualVector[ i               ] += jK_div_Mu * intg_NiCx;
            ResidualVector[ i +   mNumNodes ] += jK_div_Mu * intg_NiCy;
			ResidualVector[ i + 2*mNumNodes ] += jK_div_Mu * intg_NiCz;
        }  

		// Div element coefficient matrix ( n . E - k . E )
        Matrix<double> DMatrix( 3, 3 );

        DMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] ) - ( Next[ 0 ]*Kdir[ 0 ] );
		DMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ] ) - ( Next[ 0 ]*Kdir[ 1 ] );
		DMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ] ) - ( Next[ 0 ]*Kdir[ 2 ] );
								        
		DMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ] ) - ( Next[ 1 ]*Kdir[ 0 ] );
		DMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] ) - ( Next[ 1 ]*Kdir[ 1 ] );
		DMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ] ) - ( Next[ 1 ]*Kdir[ 2 ] );
									         
		DMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ] ) - ( Next[ 2 ]*Kdir[ 0 ] );
		DMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ] ) - ( Next[ 2 ]*Kdir[ 1 ] );
		DMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] ) - ( Next[ 2 ]*Kdir[ 2 ] );

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg_NiDx( 0.0, 0.0 );
			std::complex<double> intg_NiDy( 0.0, 0.0 );
			std::complex<double> intg_NiDz( 0.0, 0.0 );

            for( int j=0; j<mNumNodes; j++ )
			{
				intg_NiDx += ( DMatrix[ 0 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 0 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 0 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiDy += ( DMatrix[ 1 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 1 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 1 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
				intg_NiDz += ( DMatrix[ 2 ][ 0 ] * cFn[ j ][ 0 ] + DMatrix[ 2 ][ 1 ] * cFn[ j ][ 1 ] + DMatrix[ 2 ][ 2 ] * cFn[ j ][ 2 ] ) * NiNj_Matrix[ i ][ j ];
			}

            ResidualVector[ i               ] -= jK_div_Mu * intg_NiDx;
            ResidualVector[ i +   mNumNodes ] -= jK_div_Mu * intg_NiDy;
			ResidualVector[ i + 2*mNumNodes ] -= jK_div_Mu * intg_NiDz;
        }  
    }

    //**********************************************************************************************************************
    //* - Residual vector for plane waves and Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetResidualVector_Flux( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 1.0 )
        {
            GetResidualVector_Single_Wave( ResidualVector );
        }
		else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 2.0 )
		{
            GetResidualVector_Multpl_Waves( FluxValue, ResidualVector );
		}
        else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            GetResidualVector_Gauss_Beams( ResidualVector );
        }
    }

    //**********************************************************************************************************************
    //* - Stiffness matrix
    //**********************************************************************************************************************
    void GenericRobin_1st_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Setting element stiffness matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

        // Robin coefficient
        std::complex<double> cRobinCoeff( 0.0, 0.0 );

        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 )
        {
		    cRobinCoeff -= std::complex<double>( (*mProperties)( COMPLEX_IBC )[ 0 ], (*mProperties)( COMPLEX_IBC )[ 1 ] );
        }
        else
        {
            cRobinCoeff -= Calculate_jK_div_Mu();
        }
        
        // NiNj surface integrals matrix
        Matrix<double> NiNj_Matrix;

		Calculate_NiNj_Matrix( NiNj_Matrix );

        // Stiffness matrix
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                StiffMatrix[ i               ][ j               ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
            }
        }		
    }
} 