
#include "ColdPlasma.h"
#include <iostream>
#include <fstream>

namespace Kratos
{
	//***************************************************************************************************************************
	//* - Load plasma data from files.
	//***************************************************************************************************************************
	void ColdPlasma::Load_Data( String data_file_name )
	{
		// Open cold plasma properties data file in reading mode
		std::fstream dataFile( data_file_name, std::fstream::in );

		// Plasma material Id
		dataFile >> mMaterialId;

		if( dataFile.eof() )
		{
			Send_Error_Msg( "No cold plasma applied", "Please, apply plasma material to volumes." );
		}

		// Plasma geometry type:
		// 0 = Flat1D
		// 1 = Curved
		// 2 = Full3D 
		// 3 = Tensor
		dataFile >> mGeometryExtType;

		// Poloidal curvatures
		dataFile >> mPoloidalCurv_ne;
		dataFile >> mPoloidalCurv_Be;

		// First point of the 1D measurement line
		mL1_O.resize( 3 );

		dataFile >> mL1_O[ 0 ];
		dataFile >> mL1_O[ 1 ];
		dataFile >> mL1_O[ 2 ];

		// Coordinates of the center of the tokamak
		mTK_O.resize( 3 );	

		dataFile >> mTK_O[ 0 ];
		dataFile >> mTK_O[ 1 ];
		dataFile >> mTK_O[ 2 ];

		// Central axis of the tokamak
		mTK_n.resize( 3 );

		dataFile >> mTK_n[ 0 ];
		dataFile >> mTK_n[ 1 ];
		dataFile >> mTK_n[ 2 ];

		// Electron density data file
		String eDensityFile;

		dataFile >> eDensityFile;

		// External static B fields data file
		String sBfieldsFile;

		dataFile >> sBfieldsFile;

		// Densities relative to the electron density
		dataFile >> mH1_div_ne;
        dataFile >> mH2_div_ne;
        dataFile >> mH3_div_ne;

        dataFile >> mHe3_q1_div_ne;
        dataFile >> mHe3_q2_div_ne;

        dataFile >> mHe4_q1_div_ne;
        dataFile >> mHe4_q2_div_ne;

        dataFile >> mBe9_q1_div_ne;
        dataFile >> mBe9_q2_div_ne;
        dataFile >> mBe9_q3_div_ne;
        dataFile >> mBe9_q4_div_ne;

        dataFile >> mNe20_q10_div_ne; 
        dataFile >> mNi60_q20_div_ne; 
        dataFile >> mW184_q25_div_ne;

        // Matrix storage format: 
		// 0 = Full_matrix
		// 1 = Herm-Full
		// 2 = Herm-Symm 
        dataFile >> mMatrixStorageFormat;

        // Elements matrix tolerance
        dataFile >> mKij_Tol;

        // E parallel tolerance
        dataFile >> mEpar_Tol;

        // Damping mode:
		// 0 = No_damping
		// 1 = Background_ctc
		// 2 = Background_den
		// 3 = ComplexFreq_all
		// 4 = ComplexFreq_ele 
        dataFile >> mDamping_mode;

		// Damping value
		dataFile >> mDamping_value;

		// Closing cold plasma data file
        dataFile.close();

		// Complete cold plasma set-up
		if( mGeometryExtType < 2 )
		{
			Set_Tokamak_ntb();	
			Read_eDensity_File_1D( eDensityFile );
		    Read_sBfields_File_1D( sBfieldsFile );
		}
		else
		{
			Set_3D_Vector_Ne_Bxyz();
			Read_eDensity_File_3D( eDensityFile );
		    Read_sBfields_File_3D( sBfieldsFile );		
		}

		if( mGeometryExtType == 3 )
		{
		    Read_K_Tensor_File_3D( eDensityFile );
		}
	}

	//***************************************************************************************************************************
	//* - Resize 3D vectors plasma profiles.
	//***************************************************************************************************************************
	void ColdPlasma::Resize_Profile_3D_Vectors( unsigned int NewSize )
	{
	    mNe.resize( NewSize, 0.0 );
		mBx.resize( NewSize, 0.0 );
		mBy.resize( NewSize, 0.0 );
		mBz.resize( NewSize, 0.0 );

		if( mGeometryExtType == 3 ) 
		{
			int OldSize = mImp_K_Tensor.size();
			
			mImp_K_Tensor.resize( NewSize );

			for( int ki=OldSize; ki<NewSize; ki++ ) 
		    { 
			    mImp_K_Tensor[ ki ].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
		    }
		}
	}

	//***************************************************************************************************************************
	//* - Tokamak coordinates system with origen mTK_O and main axis mTK_n.
	//***************************************************************************************************************************
	void ColdPlasma::Set_Tokamak_ntb()
	{
		// Tokamak central axis unitary vector n
		double n_norm = std::sqrt( mTK_n[ 0 ]*mTK_n[ 0 ] + mTK_n[ 1 ]*mTK_n[ 1 ] + mTK_n[ 2 ]*mTK_n[ 2 ] );

		if( n_norm > 0.0 ) 
		{
			mTK_n *= ( 1.0 / n_norm );
		}
		else 
		{
			Send_Error_Msg( "Wrong tokamak vertical axis", "Please, check plasma input data." );
		}

	    // Position vector of 1D line start point from mTK_O
		Vector<double> R1D( 3 );

		R1D[ 0 ] = mL1_O[ 0 ] - mTK_O[ 0 ];
		R1D[ 1 ] = mL1_O[ 1 ] - mTK_O[ 1 ];
		R1D[ 2 ] = mL1_O[ 2 ] - mTK_O[ 2 ];

		// Unitary vector b = ( n x R1D ) / || n x R1D ||
		mTK_b.resize( 3 );

		mTK_b[ 0 ] = mTK_n[ 1 ] * R1D[ 2 ] - mTK_n[ 2 ] * R1D[ 1 ];
		mTK_b[ 1 ] = mTK_n[ 2 ] * R1D[ 0 ] - mTK_n[ 0 ] * R1D[ 2 ];
		mTK_b[ 2 ] = mTK_n[ 0 ] * R1D[ 1 ] - mTK_n[ 1 ] * R1D[ 0 ];

		double b_norm = std::sqrt( mTK_b[ 0 ]*mTK_b[ 0 ] + mTK_b[ 1 ]*mTK_b[ 1 ] + mTK_b[ 2 ]*mTK_b[ 2 ] );

		if( b_norm > 0.0 ) 
		{
			mTK_b *= ( 1.0 / b_norm );
		}
		else 
		{
			Send_Error_Msg( "Wrong L1_O_XYZ point coordinates", "Please, check plasma input data." );
		}

		// Unitary vector t = ( b x n ) / || b x n ||
		mTK_t.resize( 3 );

		mTK_t[ 0 ] = mTK_b[ 1 ] * mTK_n[ 2 ] - mTK_b[ 2 ] * mTK_n[ 1 ];
		mTK_t[ 1 ] = mTK_b[ 2 ] * mTK_n[ 0 ] - mTK_b[ 0 ] * mTK_n[ 2 ];
		mTK_t[ 2 ] = mTK_b[ 0 ] * mTK_n[ 1 ] - mTK_b[ 1 ] * mTK_n[ 0 ];

		mTK_t *= ( 1.0 / std::sqrt( mTK_t[ 0 ]*mTK_t[ 0 ] + mTK_t[ 1 ]*mTK_t[ 1 ] + mTK_t[ 2 ]*mTK_t[ 2 ] ) );

		// 1D line start point in cylindrical mTK system
		Vector<double> RZ_O( 3 );

		RZ_O[ 0 ] = R1D[ 0 ] * mTK_t[ 0 ] + R1D[ 1 ] * mTK_t[ 1 ] + R1D[ 2 ] * mTK_t[ 2 ];
		RZ_O[ 1 ] = R1D[ 0 ] * mTK_b[ 0 ] + R1D[ 1 ] * mTK_b[ 1 ] + R1D[ 2 ] * mTK_b[ 2 ];
		RZ_O[ 2 ] = R1D[ 0 ] * mTK_n[ 0 ] + R1D[ 1 ] * mTK_n[ 1 ] + R1D[ 2 ] * mTK_n[ 2 ];

		mR_start = std::sqrt( RZ_O[ 0 ]*RZ_O[ 0 ] + RZ_O[ 1 ]*RZ_O[ 1 ]);

		mZ_start = RZ_O[ 2 ];
	}

	//***************************************************************************************************************************
	//* - Transform node cartesian coordinates to tokamak coordinates.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Node_XYZ_to_TK_XYZ( Node::Pointer pNode )
	{
		// Node cartesian coordinates 
		Vector<double> OG_XYZ( 3 );

		OG_XYZ[ 0 ] = pNode->X();
		OG_XYZ[ 1 ] = pNode->Y();
		OG_XYZ[ 2 ] = pNode->Z();

		// Position vector of node from mTK_O
		Vector<double> Ptk( 3 );

		Ptk[ 0 ] = OG_XYZ[ 0 ] - mTK_O[ 0 ];
		Ptk[ 1 ] = OG_XYZ[ 1 ] - mTK_O[ 1 ];
		Ptk[ 2 ] = OG_XYZ[ 2 ] - mTK_O[ 2 ];

		// TK coordinates of the node
		Vector<double> TK_XYZ( 3 );

		TK_XYZ[ 0 ] = Ptk[ 0 ] * mTK_t[ 0 ] + Ptk[ 1 ] * mTK_t[ 1 ] + Ptk[ 2 ] * mTK_t[ 2 ];
		TK_XYZ[ 1 ] = Ptk[ 0 ] * mTK_b[ 0 ] + Ptk[ 1 ] * mTK_b[ 1 ] + Ptk[ 2 ] * mTK_b[ 2 ];
		TK_XYZ[ 2 ] = Ptk[ 0 ] * mTK_n[ 0 ] + Ptk[ 1 ] * mTK_n[ 1 ] + Ptk[ 2 ] * mTK_n[ 2 ];

		return TK_XYZ;
	}

	//***************************************************************************************************************************
	//* - Node cartesian coordinates to incR ( distance from Lo[ 0 ] ).
	//***************************************************************************************************************************
	double ColdPlasma::Node_XYZ_to_incR( Node::Pointer pNode, double PoloidalCurv )
	{
		// Node coordinates in tokamak coordinate system
		Vector<double> TK_XYZ = Node_XYZ_to_TK_XYZ( pNode );

		// Node distance to the tokamak axis
		double R;
		
		if( mGeometryExtType == 0 )
		{
			R = std::abs( TK_XYZ[ 0 ] );
		}
		else
		{
			R = std::sqrt( TK_XYZ[ 0 ]*TK_XYZ[ 0 ] + TK_XYZ[ 1 ]*TK_XYZ[ 1 ] );
		}

		// Add poloidal curvature
		if( PoloidalCurv > 0.0 )
		{
			double Z = TK_XYZ[ 2 ];

			double Ro = mR_start - PoloidalCurv;
			double Zo = mZ_start;

			double dist = std::sqrt( ( R - Ro )*( R - Ro ) + ( Z - Zo )*( Z - Zo ) );

			R = Ro + dist;
		}

        return ( mR_start - R );
	}

	//***************************************************************************************************************************
	//* - B field in XYZ cartesian coordinate system from [ modB, angB ] 1D format on pNode.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::TK_Bfield_to_XYZ( Node::Pointer pNode, double modB, double angB )
	{
		// B in ntb coordinates and xyz coordinates
		Vector<double> B_ntb( 3 );
		Vector<double> B_xyz( 3 );

		// Degrees to radians
		double pi       = 3.1415926535898;
		double angB_rad = angB * pi / 180;

		// B field in curved geometry mode
		if( mGeometryExtType == 1 )
		{
			Vector<double> N_tbn = Node_XYZ_to_TK_XYZ( pNode );

			double modN_tb = std::sqrt( N_tbn[ 0 ]*N_tbn[ 0 ] + N_tbn[ 1 ]*N_tbn[ 1 ] );

            double modB_tb = modB * std::cos( angB_rad );

			double inv_mod;

			if( modN_tb != 0.0 ) 
			{
				inv_mod = 1.0 / modN_tb;
			}
			else
			{
				inv_mod = 0.0;
			}

			B_ntb[ 0 ] = + modB_tb * N_tbn[ 1 ] * inv_mod;
			B_ntb[ 1 ] = - modB_tb * N_tbn[ 0 ] * inv_mod;
			B_ntb[ 2 ] = + modB    * std::sin( angB_rad );
		}
		// B field in flat geometry mode
		else
		{
			B_ntb[ 0 ] = + 0.0;
			B_ntb[ 1 ] = - modB * std::cos( angB_rad );
			B_ntb[ 2 ] = + modB * std::sin( angB_rad );
		}

		// Add poloidal curvature to Be field
		if( mPoloidalCurv_Be > 0.0 ) 
		{
			Add_PoloidalCurv_to_TKBfield( pNode, B_ntb );
		}

		// B in ntb CS to xyz CS
		B_xyz[ 0 ] = B_ntb[ 0 ] * mTK_t[ 0 ] + B_ntb[ 1 ] * mTK_b[ 0 ] + B_ntb[ 2 ] * mTK_n[ 0 ];
		B_xyz[ 1 ] = B_ntb[ 0 ] * mTK_t[ 1 ] + B_ntb[ 1 ] * mTK_b[ 1 ] + B_ntb[ 2 ] * mTK_n[ 1 ];
		B_xyz[ 2 ] = B_ntb[ 0 ] * mTK_t[ 2 ] + B_ntb[ 1 ] * mTK_b[ 2 ] + B_ntb[ 2 ] * mTK_n[ 2 ];

		return B_xyz;
	}

	//***************************************************************************************************************************
	//* - Add poloidal curvature to B ext.
	//***************************************************************************************************************************
	void ColdPlasma::Add_PoloidalCurv_to_TKBfield( Node::Pointer pNode, Vector<double>& B_ntb )
	{
		// Node coordinates in TK system
		Vector<double> rNode = Node_XYZ_to_TK_XYZ( pNode );

		// Curvature center coordinates in TK sytem
		Vector<double> rCurC( 3 );

		// Curved geometry mode
		if( mGeometryExtType == 1 )
		{
		    double modN_tb = std::sqrt( rNode[ 0 ]*rNode[ 0 ] + rNode[ 1 ]*rNode[ 1 ] );

		    double inv_modN_tb;

			if( modN_tb != 0.0 ) 
			{
				inv_modN_tb = 1.0 / modN_tb;
			}
			else  
			{
				inv_modN_tb = 0.0;
			}

			rCurC[ 0 ] = ( mR_start - mPoloidalCurv_Be ) * rNode[ 0 ] * inv_modN_tb;
			rCurC[ 1 ] = ( mR_start - mPoloidalCurv_Be ) * rNode[ 1 ] * inv_modN_tb;
			rCurC[ 2 ] = ( mZ_start );
		}
		// Flat geometry mode
		else
		{
			rNode[ 1 ] = 0.0;

			rCurC[ 0 ] = mR_start - mPoloidalCurv_Be;
			rCurC[ 1 ] = 0.0;
			rCurC[ 2 ] = mZ_start;
		}

		// Node position vector from curvature center
		Vector<double> rN_CC( 3 );

		rN_CC[ 0 ] = rNode[ 0 ] - rCurC[ 0 ];
		rN_CC[ 1 ] = rNode[ 1 ] - rCurC[ 1 ];
		rN_CC[ 2 ] = rNode[ 2 ] - rCurC[ 2 ];

		// Vector perpendicular to rN_CC and mTK_n axis ( P_rN_n = rN_CC x mTK_n ) 
		Vector<double> P_rN_n( 3 );

		P_rN_n[ 0 ] = + rN_CC[ 1 ];
		P_rN_n[ 1 ] = - rN_CC[ 0 ];
		P_rN_n[ 2 ] = 0.0;

		// Vector perpendicular to P_rN_n and rN_CC ( P_rN_CC = P_rN_n x rN_CC )
		Vector<double> P_rN_CC( 3 );

		P_rN_CC[ 0 ] = P_rN_n[ 1 ] * rN_CC[ 2 ] - P_rN_n[ 2 ] * rN_CC[ 1 ];
		P_rN_CC[ 1 ] = P_rN_n[ 2 ] * rN_CC[ 0 ] - P_rN_n[ 0 ] * rN_CC[ 2 ];
		P_rN_CC[ 2 ] = P_rN_n[ 0 ] * rN_CC[ 1 ] - P_rN_n[ 1 ] * rN_CC[ 0 ];

		// Change direction after crossing mTK_n axis
		if( rN_CC[ 0 ] < 0 )
		{
		    P_rN_CC *= -1;
		}

		// Change direction after crossing mTK_b axis
		if( rCurC[ 0 ] < 0 )
		{
		    P_rN_CC *= -1;
		}

		// Redefine P_rN_CC if rN_CC is parallel to mTK_n
		if( ( rN_CC[ 0 ] == 0.0 ) && ( rN_CC[ 1 ] == 0.0 ) )
		{
		    P_rN_CC = rCurC;
			
			if( rN_CC[ 2 ] > 0.0 ) 
			{
				P_rN_CC *= -1;
			}
		}
		
		// P_rN_CC vector normalization
		double modPrNCC = std::sqrt( P_rN_CC[ 0 ]*P_rN_CC[ 0 ] + P_rN_CC[ 1 ]*P_rN_CC[ 1 ] + P_rN_CC[ 2 ]*P_rN_CC[ 2 ] );

	    if( modPrNCC > 0.0 ) 
		{
			P_rN_CC *= ( 1.0 / modPrNCC );
		}
		
		// Projection of Bz on the normalized vector P_rN_CC
		B_ntb[ 0 ] += B_ntb[ 2 ] * P_rN_CC[ 0 ];
		B_ntb[ 1 ] += B_ntb[ 2 ] * P_rN_CC[ 1 ];
		B_ntb[ 2 ]  = B_ntb[ 2 ] * P_rN_CC[ 2 ];
	}

	//***************************************************************************************************************************
	//* - Lower and upper indexes in mLo containing incR:
	//*   index_in_Lo[ 0 ] = lower_index
	//*   index_in_Lo[ 1 ] = upper_index
	//***************************************************************************************************************************
	Vector<int> ColdPlasma::Get_Indexes_In_Lo( double incR )
	{
	    // Closest value of incR in mLo
		Vector<double>::iterator it_mLo;

		// Upper index ( it_mLo - mLo.begin() )
		it_mLo = std::lower_bound( mLo.begin(), mLo.end(), incR );

		// Indexes containing incR between them
		Vector<int> index_in_Lo( 2 );

		if( it_mLo == mLo.begin() )
		{
			index_in_Lo[ 0 ] = 0;
			index_in_Lo[ 1 ] = 0;
		}
		else if( it_mLo == mLo.end() )
		{
			index_in_Lo[ 0 ] = mLo.size() - 1;
			index_in_Lo[ 1 ] = mLo.size() - 1;
		}
		else
		{
			index_in_Lo[ 0 ] = it_mLo - mLo.begin() - 1;
			index_in_Lo[ 1 ] = it_mLo - mLo.begin();
		}

		return index_in_Lo;
	}

	//***************************************************************************************************************************
	//* - Bext in XYZ cartesian components from vectors mB_mod and mB_ang by linear interpolation.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Get_Bext_1D( double incR, Vector<int>& index, Node::Pointer pNode )
	{
		double modB, angB;

		if( index[ 0 ] == index[ 1 ] )
		{
			modB = mB_mod[ index[ 0 ] ];
			angB = mB_ang[ index[ 0 ] ];
		}
		else
		{
			double alpha = ( incR - mLo[ index[ 0 ] ] ) / std::abs( mLo[ index[ 1 ] ] - mLo[ index[ 0 ] ] );

			modB = alpha * mB_mod[ index[ 1 ] ] + ( 1.0 - alpha ) * mB_mod[ index[ 0 ] ];
			angB = alpha * mB_ang[ index[ 1 ] ] + ( 1.0 - alpha ) * mB_ang[ index[ 0 ] ];
		}

		return TK_Bfield_to_XYZ( pNode, modB, angB );
	}

	//****************************************************************************************************************************
	//* - Electron density from vector mNe by linear interpolation.
	//****************************************************************************************************************************
	double ColdPlasma::Get_eDensity_1D( double incR, Vector<int>& index )
	{
		double eDensityOnNode;

		if( index[ 0 ] == index[ 1 ] )
		{
			eDensityOnNode = mNe[ index[ 0 ] ];
		}
		else
		{
			double alpha = ( incR - mLo[ index[ 0 ] ] ) / std::abs( mLo[ index[ 1 ] ] - mLo[ index[ 0 ] ] );

			eDensityOnNode = alpha * mNe[ index[ 1 ] ] + ( 1.0 - alpha ) * mNe[ index[ 0 ] ];
		}

		return eDensityOnNode;
	}

	//***************************************************************************************************************************
	//* - eDensity on pNode.
	//***************************************************************************************************************************
	double ColdPlasma::Get_eDensity_OnNode( Node::Pointer pNode )
	{
		double eDensity = 0.0;

	    if( mGeometryExtType < 2 )
	    {
		    double      incRe = Node_XYZ_to_incR ( pNode, mPoloidalCurv_ne );
		    Vector<int> indxe = Get_Indexes_In_Lo( incRe );
		    
		    eDensity = Get_eDensity_1D( incRe, indxe );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( indxp < mNe.size() )
			{
                eDensity = mNe[ indxp ];
			}			
		}

		return eDensity;
	}

	//***************************************************************************************************************************
	//* - Bext on pNode.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Get_Bext_OnNode( Node::Pointer pNode )
	{
		Vector<double> Bext( 3, 0.0 );

	    if( mGeometryExtType < 2 )
	    {
		    double      incRb = Node_XYZ_to_incR ( pNode, mPoloidalCurv_Be );
		    Vector<int> indxb = Get_Indexes_In_Lo( incRb );
		    
		    Bext = Get_Bext_1D( incRb, indxb, pNode );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( indxp < mBx.size() )
			{
                Bext[ 0 ] = mBx[ indxp ];
			    Bext[ 1 ] = mBy[ indxp ];
			    Bext[ 2 ] = mBz[ indxp ];
			}			
		}

		return Bext;
	}

    //***************************************************************************************************************************
    // - eDensity and Bext on pNode.
    //***************************************************************************************************************************
    void ColdPlasma::Get_eDensity_Bext_OnNode( double& eDensity, Vector<double>& Bext, Node::Pointer pNode )
    {
	    eDensity = 0.0;
			
		Bext.resize( 3, 0.0 );	    
		
		if( mGeometryExtType < 2 )
	    {
		    double      incRe = Node_XYZ_to_incR ( pNode, mPoloidalCurv_ne );
		    Vector<int> indxe = Get_Indexes_In_Lo( incRe );

            eDensity = Get_eDensity_1D( incRe, indxe );
		    
		    double      incRb = Node_XYZ_to_incR ( pNode, mPoloidalCurv_Be );
		    Vector<int> indxb = Get_Indexes_In_Lo( incRb );

		    Bext = Get_Bext_1D( incRb, indxb, pNode );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( indxp < mNe.size() )
			{
                eDensity = mNe[ indxp ];
			}

			if( indxp < mBx.size() )
			{
                Bext[ 0 ] = mBx[ indxp ];
			    Bext[ 1 ] = mBy[ indxp ];
			    Bext[ 2 ] = mBz[ indxp ];
			}			
		}
    }

    //***************************************************************************************************************************
	// - eDensity, Bext and [ S, D, P, R, L ] components of the permittivity tensor on pNode.
    //***************************************************************************************************************************
	void ColdPlasma::Get_AllPlasmaParameters_OnNode( double& eDensity, Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq )
	{
		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
	}

    //***************************************************************************************************************************
    // - Bext and [ S, D, P, R, L ] components of the permittivity tensor on pNode.
    //***************************************************************************************************************************
    void ColdPlasma::Get_Bext_SDPRL_OnNode( Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq )
    {
        double eDensity; 

		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
    }

	//***************************************************************************************************************************
	//* - [ S, D, P, R, L ] components of the permittivity tensor on pNode.
	//***************************************************************************************************************************
    ColdPlasma::ComplexVector ColdPlasma::Get_SDPRL_OnNode( Node::Pointer pNode, double wFreq )
	{
		double eDensity;

		Vector<double> Bext;

		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        return Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
	}

	//***************************************************************************************************************************
	//* - [ S, D, P, R, L ] components of the permittivity tensor.
	//***************************************************************************************************************************
    ColdPlasma::ComplexVector ColdPlasma::Calculate_PlasmaPermittivity_SDPRL( double eDensity, Vector<double>& Bext, double wFreq )
	{
		// Module of Bext
		double modBext  = std::sqrt( Bext[ 0 ]*Bext[ 0 ] + Bext[ 1 ]*Bext[ 1 ] + Bext[ 2 ]*Bext[ 2 ] );

		// Useful constants
		double eo = 8.8542e-12; // Vacuum electric permittivity
		double Qe = 1.6022e-19; // Electron charge
		double Me = 9.1095e-31; // Electron mass
        double Wf = wFreq;      // Problem frequency (angular)

        // Complex unit
        std::complex<double> cUnit( 0.0, 1.0 );

        // Atomic Mass Unit in Kg
        double uKg   = 1.660540199e-27;

        double MH1   = 1.0078250 * uKg; // Hidrogen mass
        double MH2   = 2.0141018 * uKg; // Deuterium mass
        double MH3   = 3.0160493 * uKg; // Tritium mass
        double MHe3  = 3.0160293 * uKg; // Helium 3 mass
        double MHe4  = 4.0026032 * uKg; // Helium 4 mass
        double MBe9  = 9.0121822 * uKg; // Beryllium 9 mass
        double MNe20 = 19.992440 * uKg; // Neon 20 mass
        double MNi60 = 59.930786 * uKg; // Nickel 60 mass
        double MW184 = 183.95093 * uKg; // Tungsten 184 mass

        // Total number of species
        int numSpecies = 15;

		// Defining species properties (charge, mass and density)
		Vector<double> Qs( numSpecies );
		Vector<double> Ms( numSpecies );
		Vector<double> Ns( numSpecies );

		// Electron
		Qs[ 0 ] =-Qe;
		Ms[ 0 ] = Me;
		Ns[ 0 ] = eDensity;

		// Hydrogen (H1)
		Qs[ 1 ] = Qe;
		Ms[ 1 ] = MH1;
		Ns[ 1 ] = eDensity * mH1_div_ne;

        // Deuterium (H2)
		Qs[ 2 ] = Qe;
		Ms[ 2 ] = MH2;
		Ns[ 2 ] = eDensity * mH2_div_ne;

        // Tritium (H3)
		Qs[ 3 ] = Qe;
		Ms[ 3 ] = MH3;
		Ns[ 3 ] = eDensity * mH3_div_ne;

        // Helium 3 (He3), Q=+1
        Qs[ 4 ] = Qe;
        Ms[ 4 ] = MHe3;
        Ns[ 4 ] = eDensity * mHe3_q1_div_ne;

        // Helium 3 (He3), Q=+2
        Qs[ 5 ] = 2.0 * Qe;
        Ms[ 5 ] = MHe3;
        Ns[ 5 ] = eDensity * mHe3_q2_div_ne;

        // Helium 4 (He4), Q=+1
        Qs[ 6 ] = Qe;
        Ms[ 6 ] = MHe4;
        Ns[ 6 ] = eDensity * mHe4_q1_div_ne;

        // Helium 4 (He4), Q=+2
        Qs[ 7 ] = 2.0 * Qe;
        Ms[ 7 ] = MHe4;
        Ns[ 7 ] = eDensity * mHe4_q2_div_ne;

        // Beryllium 9 (Be9), Q=+1
        Qs[ 8 ] = Qe;
        Ms[ 8 ] = MBe9;
        Ns[ 8 ] = eDensity * mBe9_q1_div_ne;

        // Beryllium 9 (Be9), Q=+2
        Qs[ 9 ] = 2.0 * Qe;
        Ms[ 9 ] = MBe9;
        Ns[ 9 ] = eDensity * mBe9_q2_div_ne;

        // Beryllium 9 (Be9), Q=+3
        Qs[ 10 ] = 3.0 * Qe;
        Ms[ 10 ] = MBe9;
        Ns[ 10 ] = eDensity * mBe9_q3_div_ne;

        // Beryllium 9 (Be9), Q=+4
        Qs[ 11 ] = 4.0 * Qe;
        Ms[ 11 ] = MBe9;
        Ns[ 11 ] = eDensity * mBe9_q4_div_ne;

       // Neon 10 (Ne10), Q=+10
        Qs[ 12 ] = 10.0 * Qe;
        Ms[ 12 ] = MNe20;
        Ns[ 12 ] = eDensity * mNe20_q10_div_ne;

        // Nickel 60 (Ni60), Q=+20
        Qs[ 13 ] = 20.0 * Qe;
        Ms[ 13 ] = MNi60;
        Ns[ 13 ] = eDensity * mNi60_q20_div_ne;

        // Tungsten 184 (W184), Q=+25
        Qs[ 14 ] = 25.0 * Qe;
        Ms[ 14 ] = MW184;
        Ns[ 14 ] = eDensity * mW184_q25_div_ne;

        // Plasma frequencies
		Vector<double> W2_ps( numSpecies );

		for( int i=0; i<numSpecies; i++ )
		{
			W2_ps[ i ] = Qs[ i ] * Qs[ i ] * Ns[ i ] / ( eo * Ms[ i ] );
		}

		// Ciclotron frequencies
		Vector<double> Omega_s( numSpecies );

		for( int i=0; i<numSpecies; i++ )
		{
			Omega_s[ i ] = Qs[ i ] * modBext / Ms[ i ];
		}

		// Collisional frequencies 
		Vector< std::complex<double> > W_col( numSpecies );

		if( mDamping_mode == 3 )
		{
		    for( int i=0; i<numSpecies; i++ )
		    {
			    W_col[ i ] = std::complex<double>( Wf, mDamping_value );
		    }		
		}
		else if( mDamping_mode == 4 )
		{
		    W_col[ 0 ] = std::complex<double>( Wf, mDamping_value );
			
			for( int i=1; i<numSpecies; i++ )
		    {
			    W_col[ i ] = std::complex<double>( Wf, 0.0 );
		    }		
		}
		else
		{
		    for( int i=0; i<numSpecies; i++ )
		    {
			    W_col[ i ] = std::complex<double>( Wf, 0.0 );
		    }		
		}

        // Computing R
        std::complex<double> R( 1.0, 0.0 );

        for( int i=0; i<numSpecies; i++ )
		{
			if( std::abs( Wf + Omega_s[ i ] ) > 0.0 )
			{
				R -= W2_ps[ i ] / ( Wf * ( W_col[ i ] + Omega_s[ i ] ) );
			}
			else
			{
				Send_Error_Msg( "R plasma parameter infinite", "Please, please check input data." );
			}
		}

        // Computing L
        std::complex<double> L( 1.0, 0.0 );

        for( int i=0; i<numSpecies; i++ )
		{
			if( std::abs( Wf - Omega_s[ i ] ) > 0.0 )
			{
				L -= W2_ps[ i ] / ( Wf * ( W_col[ i ] - Omega_s[ i ] ) );
			}
			else
			{
				Send_Error_Msg( "L plasma parameter infinite", "Please, please check input data." );
			}
		}

        // Computing P
        std::complex<double> P( 1.0, 0.0 );

        for( int i=0; i<numSpecies; i++ )
		{
			P -= W2_ps[ i ] / ( Wf * W_col[ i ] );
		}

		// Adding damping in 'Background_const' damping mode
		if( mDamping_mode == 1 )
		{
		    R += cUnit * mDamping_value;
		    L += cUnit * mDamping_value;
		    P += cUnit * mDamping_value;
		}

		// Adding damping in 'Background_denst' damping mode
		if( mDamping_mode == 2 )
		{
			if( eDensity > 0.0 )
			{
		        R += cUnit * mDamping_value;
		        L += cUnit * mDamping_value;
		        P += cUnit * mDamping_value;
			}
		}

		// Components S, D and P of the plasma permittivity tensor (Stix notation)
		// Components R, L and P from the diagonal plasma permittivity tensor
        ComplexVector SDPRL( 5 );

		SDPRL[ 0 ] = 0.5 * ( R + L );
	    SDPRL[ 1 ] = 0.5 * ( R - L );
	    SDPRL[ 2 ] = P;
		SDPRL[ 3 ] = R;
		SDPRL[ 4 ] = L;

        return SDPRL;
	}

    //***************************************************************************************************************************
    //* - Rotation matrix.
    //***************************************************************************************************************************
    void ColdPlasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz )
    {
         // Unit vectors
         Vector<double> n( 3 ), t( 3 ), b( 3 );

         // b = B / |B|
         double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

         b[ 0 ] = Bx / Bnorm;
         b[ 1 ] = By / Bnorm;
         b[ 2 ] = Bz / Bnorm;

         // n _|_ b
         double nnorm = std::sqrt( b[ 0 ]*b[ 0 ] + b[ 1 ]*b[ 1 ] );

         if( nnorm > 0.0 )
         {
             n[ 0 ] = b[ 1 ] / nnorm;
             n[ 1 ] =-b[ 0 ] / nnorm;
             n[ 2 ] = 0.0;
         }
         else
         {
             n[ 0 ] = b[ 2 ];
             n[ 1 ] = 0.0;
             n[ 2 ] = 0.0;
         }

         // t = b x n
         t[ 0 ] = b[ 1 ]*n[ 2 ] - b[ 2 ]*n[ 1 ];
         t[ 1 ] = b[ 2 ]*n[ 0 ] - b[ 0 ]*n[ 2 ];
         t[ 2 ] = b[ 0 ]*n[ 1 ] - b[ 1 ]*n[ 0 ];

         // Rotation matrix
         // T = [ nx tx bx
         //       ny ty by
         //       nz tz bz ]
         T.Resize( 3, 3 );

         for( int row=0; row<3; row++ )
         {
             T[ row ][ 0 ] = std::complex<double>( n[ row ] );
             T[ row ][ 1 ] = std::complex<double>( t[ row ] );
             T[ row ][ 2 ] = std::complex<double>( b[ row ] );
         }
    }

	//***************************************************************************************************************************
	// - Permittivity tensor on pNode.
	//***************************************************************************************************************************
	void ColdPlasma::Get_PermittivityTensor_OnNode( Matrix<std::complex<double> >& TEn, Node::Pointer pNode, double wFreq )
	{
		// Permittivity tensor in a node
		TEn.Resize( 3, 3 );

		// Relative permittivity tensor in b = B/|B| coordinates
		Matrix<std::complex<double> > K( 3, 3 );

		// Bext on node
		Vector<double> Bext = Get_Bext_OnNode( pNode );

		// Constants
        double eo = 8.8541878176e-12;

		// Get K relative in b = B/|B| coordinates
		if( mGeometryExtType == 3 )
		{
		    K = mImp_K_Tensor[ pNode->Id() - 1 ];
		}
		else
		{
		    double eDensity = Get_eDensity_OnNode( pNode );

		    ColdPlasma::ComplexVector SDPRL;

            SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
			
			std::complex<double> S = SDPRL[ 0 ];
			std::complex<double> D = SDPRL[ 1 ];
			std::complex<double> P = SDPRL[ 2 ];

			std::complex<double> cZero( 0.0, 0.0 );
			std::complex<double> cUnit( 0.0, 1.0 );

		    K[ 0 ][ 0 ] = S        ; K[ 0 ][ 1 ] =-cUnit * D; K[ 0 ][ 2 ] = cZero;
			K[ 1 ][ 0 ] = cUnit * D; K[ 1 ][ 1 ] = S        ; K[ 1 ][ 2 ] = cZero;
			K[ 2 ][ 0 ] = cZero    ; K[ 2 ][ 1 ] = cZero    ; K[ 2 ][ 2 ] = P    ;
	    }

		// K absolute
		K *= eo;

		// Rotation matrix
		Matrix< std::complex<double> > T;

		Rotation_Matrix( T, Bext[0], Bext[1], Bext[2] );

		// TEn = T * K * Tt
		// Be carefull!!, T.Transpose() changes T permanently
		TEn = T   * K;
		TEn = TEn * T.Transpose();
	}

    //***************************************************************************************************************************
	// - Conductivity tensor on pNode.
	//***************************************************************************************************************************
	void ColdPlasma::Get_ConductivityTensor_OnNode( Matrix<std::complex<double> >& TCn, Node::Pointer pNode, double wFreq )
	{
		// Coductivity tensor in a node
		TCn.Resize( 3, 3 );

        // Permittivity tensor in a node
        Matrix< std::complex<double> > TEn;

        Get_PermittivityTensor_OnNode( TEn, pNode, wFreq );

        // Constants
        double eo = 8.8541878176e-12;

        std::complex<double> iw( 0.0, wFreq );

        TCn[ 0 ][ 0 ] = -iw * ( TEn[ 0 ][ 0 ] - eo );
		TCn[ 1 ][ 0 ] = -iw * ( TEn[ 1 ][ 0 ]      );
		TCn[ 2 ][ 0 ] = -iw * ( TEn[ 2 ][ 0 ]      );

        TCn[ 0 ][ 1 ] = -iw * ( TEn[ 0 ][ 1 ]      );
		TCn[ 1 ][ 1 ] = -iw * ( TEn[ 1 ][ 1 ] - eo );
		TCn[ 2 ][ 1 ] = -iw * ( TEn[ 2 ][ 1 ]      );

        TCn[ 0 ][ 2 ] = -iw * ( TEn[ 0 ][ 2 ]      );
		TCn[ 1 ][ 2 ] = -iw * ( TEn[ 1 ][ 2 ]      );
		TCn[ 2 ][ 2 ] = -iw * ( TEn[ 2 ][ 2 ] - eo );
	}

	//***************************************************************************************************************************
	//* - Read electron density file in 1D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_eDensity_File_1D( String FileName )
	{
		std::fstream eDensityFile( FileName, std::fstream::in );

		if( !eDensityFile.is_open() )
		{
			Send_Error_Msg( "Error opening 1D electron density file", "Please, check file path." );
		}
		
        String FileLine;

		double Lo, Ne;
		
        while( std::getline( eDensityFile, FileLine ) ) 
		{ 
		    if( isdigit( FileLine[ 0 ] ) )
			{
			    std::stringstream LineStrm( FileLine );
				
			    LineStrm >> Lo; 
				LineStrm >> Ne;
				
                mLo.push_back( Lo );
                mNe.push_back( Ne );	
		    }
        } 

		eDensityFile.close();
	}

	//***************************************************************************************************************************
	//* - Set 3D Ne and Bxyz vectors size.
	//***************************************************************************************************************************
	void ColdPlasma::Set_3D_Vector_Ne_Bxyz()
	{
        int MaxNodeId = 0;
           
        std::fstream Nodes_File( mBase_FileName + "-1.dat", std::fstream::in );  

        // Line from volume element file
		std::string line;

		// Node Id and node coordinates
        int NodeId; double X, Y, Z;

        while( std::getline( Nodes_File, line ) )
        {
            // Ignore comments
            if( line.substr( 0, 2 ) == "//" ) continue;

            // Read formated line
            std::sscanf( line.c_str(), "No[%d] = p(%lf,%lf,%lf);", &NodeId, &X, &Y, &Z );

			// Update MaxNodeId
            if( NodeId >= MaxNodeId ) MaxNodeId = NodeId;
        }

        Nodes_File.close();

		Resize_Profile_3D_Vectors( MaxNodeId );
	}

	//***************************************************************************************************************************
	//* - Read electron density file in 3D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_eDensity_File_3D( String FileName )
	{
		std::fstream eDensityFile( FileName, std::fstream::in );

		if( !eDensityFile.is_open() )
		{
			Send_Error_Msg( "Error opening 3D electron density file", "Please, check file path." );
		}
		
		String FileLine;

		int nId; double nNe;

		while( std::getline( eDensityFile, FileLine ) ) 
		{ 
		    if( isdigit( FileLine[ 0 ] ) )
			{
			    std::stringstream LineStrm( FileLine );
				
			    LineStrm >> nId; 
				LineStrm >> nNe;

				if( nId > mNe.size() )
				{
					Send_Error_Msg( "Index out of range in density file", "Please, check density file." );
				}
				
                mNe[ nId - 1 ] = nNe;	
		    }
        } 

		eDensityFile.close();
	}

	//***************************************************************************************************************************
	//* - Read external static B fields file in 1D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_sBfields_File_1D( String FileName )
	{
		std::fstream sBfieldsFile( FileName, std::fstream::in );

		if( !sBfieldsFile.is_open() )
		{
			Send_Error_Msg( "Error opening 1D external static B fields file", "Please, check file path." );
		}

        String FileLine;

		double Lo, modB, angB;
		
        while( std::getline( sBfieldsFile, FileLine ) )  
		{ 
		    if( isdigit( FileLine[ 0 ] ) )   
			{
			    std::stringstream LineStrm( FileLine );
				
			    LineStrm >> Lo  ;
			    LineStrm >> modB;
				LineStrm >> angB;
				
                mB_mod.push_back( modB );
                mB_ang.push_back( angB );	
		    }
        } 

		sBfieldsFile.close();
	}

	//***************************************************************************************************************************
	//* - Read external static B fields file in 3D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_sBfields_File_3D( String FileName )
	{
		std::fstream sBfieldsFile( FileName, std::fstream::in );

		if( !sBfieldsFile.is_open() )
		{
			Send_Error_Msg( "Error opening 3D external static B fields file", "Please, check file path." );
		}

		String FileLine;

		int nId; double nBx, nBy, nBz;

        while( std::getline( sBfieldsFile, FileLine ) )  
		{ 
		    if( isdigit( FileLine[ 0 ] ) )   
			{
			    std::stringstream LineStrm( FileLine );
				
			    LineStrm >> nId;
			    LineStrm >> nBx;
				LineStrm >> nBy;
				LineStrm >> nBz;

				if( nId > mBx.size() )
				{
					Send_Error_Msg( "Index out of range in Be field file", "Please, check Be field file." );
				}

				mBx[ nId - 1 ] = nBx;
                mBy[ nId - 1 ] = nBy;
				mBz[ nId - 1 ] = nBz;	
		    }
        } 

		sBfieldsFile.close();
	}

	//***************************************************************************************************************************
	//* - Read imported K tensor from file.
	//***************************************************************************************************************************
	void ColdPlasma::Read_K_Tensor_File_3D( String FileName )
	{
		// Imported K tensor file
		std::ifstream K_Tensor_file( "Tensor_K_CPGE.bin", std::ios::binary );

		if( !K_Tensor_file.is_open() )
		{
			Send_Error_Msg( "Error opening imported K tensor file", "Please, check file location." );
		}

		// Read ne density file 3D for node index data 
		std::fstream Ne_file( FileName, std::fstream::in );

		// Calculate binary matrix size
        int VectorSize         = 3; 
        int MatrixSize         = VectorSize * VectorSize; 
        int Chunck_Size_Matrix = MatrixSize * sizeof( std::complex<double> );

		// Resize K tensors vector
		mImp_K_Tensor.resize( mNe.size() );

		for( int ki=0; ki<mNe.size(); ki++ ) 
		{ 
			mImp_K_Tensor[ ki ].Resize( VectorSize, VectorSize, std::complex<double>( 0.0, 0.0 ) );
		}

		// Read imported K tensor from file
		std::string line;	

		int nId; double nNe;

        while( std::getline( Ne_file, line ) ) 
        {
			std::stringstream LineStrm( line );
				
			LineStrm >> nId; 
			LineStrm >> nNe;

            Matrix< std::complex<double> > K_Tensor;
			
			K_Tensor.Resize( VectorSize, VectorSize, std::complex<double>(0.0, 0.0) );

            Vector< std::complex<double> > Kb( MatrixSize, std::complex<double>(0.0, 0.0) );

            K_Tensor_file.read( reinterpret_cast<char*>( Kb.data() ), Chunck_Size_Matrix );

            if( K_Tensor_file.gcount() == Chunck_Size_Matrix )
            {
                for( int i=0; i<VectorSize; i++ )
                {
                    for( int j=0; j<VectorSize; j++ )
                    {
                        K_Tensor[ i ][ j ] = Kb[ i * VectorSize + j ];
                    }
                }
            }

			mImp_K_Tensor[ nId - 1 ] = K_Tensor;

            Kb      .FreeData();
            K_Tensor.FreeData();
        }

        Ne_file      .close();
        K_Tensor_file.close();
	}

	//***************************************************************************************************************************
	//* - Linear interpolation of plasma profiles to high-order nodes.
	//***************************************************************************************************************************
	void ColdPlasma::Linear_Profile_Interpolate( Vector<int> HONodesId )
	{
		// Interpolate plasma profile to 2nd order middle nodes
		if( HONodesId.size() == 10 )
		{
		    // 2nd order edge node Ids
		    int iniNId, midNId, endNId; 
		    
		    // Edge[ 0 , 1 ]
		    iniNId = HONodesId[ 0 ] - 1; 
		    midNId = HONodesId[ 4 ] - 1; 
		    endNId = HONodesId[ 1 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}

		    // Edge[ 1 , 2 ]
		    iniNId = HONodesId[ 1 ] - 1; 
		    midNId = HONodesId[ 5 ] - 1; 
		    endNId = HONodesId[ 2 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}

		    // Edge[ 2 , 0 ]
		    iniNId = HONodesId[ 2 ] - 1; 
		    midNId = HONodesId[ 6 ] - 1; 
		    endNId = HONodesId[ 0 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 0 , 3 ]
		    iniNId = HONodesId[ 0 ] - 1; 
		    midNId = HONodesId[ 7 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 1 , 3 ]
		    iniNId = HONodesId[ 1 ] - 1; 
		    midNId = HONodesId[ 8 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 2 , 3 ]
		    iniNId = HONodesId[ 2 ] - 1; 
		    midNId = HONodesId[ 9 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		}
		// Interpolate plasma profile to 3sb surface middle nodes
		else if ( HONodesId.size() == 17 )
		{
		    // 3sb surface vertices nodes Ids
		    int vNId0, vNId1, vNId2; 

			// 3sb surface inner nodes Ids
			int sNId0, sNId1, sNId2;  
		    
		    // Face[ 0 , 1 , 2 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 1 ] - 1; 
		    vNId2 = HONodesId[ 2 ] - 1;

			sNId0 = HONodesId[ 4 ] - 1; 
		    sNId1 = HONodesId[ 5 ] - 1; 
		    sNId2 = HONodesId[ 6 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 0 , 1 , 3 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 1 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 7 ] - 1; 
		    sNId1 = HONodesId[ 8 ] - 1; 
		    sNId2 = HONodesId[ 9 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 1 , 2 , 3 ]
		    vNId0 = HONodesId[ 1 ] - 1; 
		    vNId1 = HONodesId[ 2 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 10 ] - 1; 
		    sNId1 = HONodesId[ 11 ] - 1; 
		    sNId2 = HONodesId[ 12 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 0 , 2 , 3 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 2 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 13 ] - 1; 
		    sNId1 = HONodesId[ 14 ] - 1; 
		    sNId2 = HONodesId[ 15 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}
		}
	}

	//***************************************************************************************************************************
	//* - Print error messages and exit ERMES.
	//***************************************************************************************************************************
	void ColdPlasma::Send_Error_Msg( String ErrorMsg, String CheckMsg )
	{
		std::cout << std::endl;

		std::cout << "!!!!!!!!!! " << ErrorMsg << " !!!!!!!!!!" << std::endl;

		std::cout << std::endl;

		std::cout << CheckMsg << std::endl;

		std::cout << std::endl;

		std::cout << "ERMES analysis finished." << std::endl;

		exit( EXIT_SUCCESS );
	}

	//***************************************************************************************************************************
	//* - Clear cold plasma properties.
	//***************************************************************************************************************************
	ColdPlasma::~ColdPlasma()
	{
	    mLo.clear(); Vector<double>().swap( mLo );
		mNe.clear(); Vector<double>().swap( mNe );
	    
	    mBx.clear(); Vector<double>().swap( mBx );
	    mBy.clear(); Vector<double>().swap( mBy );
	    mBz.clear(); Vector<double>().swap( mBz );

		mL1_O.clear(); Vector<double>().swap( mL1_O );
	    mTK_O.clear(); Vector<double>().swap( mTK_O );	

		mTK_n.clear(); Vector<double>().swap( mTK_n );
	    mTK_t.clear(); Vector<double>().swap( mTK_t );
	    mTK_b.clear(); Vector<double>().swap( mTK_b );

	    mB_mod.clear(); Vector<double>().swap( mB_mod );
	    mB_ang.clear(); Vector<double>().swap( mB_ang );

		mImp_K_Tensor.clear(); Vector< Matrix<std::complex<double> > >().swap( mImp_K_Tensor ); 
	}
} 
