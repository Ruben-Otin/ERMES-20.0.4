#include <math.h>
#include <numeric>
#include <complex>
#include <omp.h>

#include "complex_matrix_operations.h"

namespace Kratos
{
    ////*********************************************************************************
    ////*  Y = A*X with A symmetric matrix
    ////*********************************************************************************
    //void cMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    //{
    //    ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
    //    ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();
    //
    //    ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
    //    ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;
    //
    //    ComplexVectorType::iterator y_iterator = Y.begin();
    //    ComplexVectorType::iterator x_iterator = X.begin();
    //
    //    std::fill(Y.begin(), Y.end(), ComplexType());
    //    
    //    ComplexType aij;
    //
    //    int i_row = 0;
    //    int i_col = 0;
    //
    //    while(row_iterator != end_iterator)
    //    {
    //        data_iterator     = row_iterator->begin();
    //        end_data_iterator = row_iterator->end();
    //
    //        ComplexType result(0.00,0.00);
    //
    //        while(data_iterator != end_data_iterator)
    //        {
    //            i_col = data_iterator->first;
    //            aij   = data_iterator->second;
    //
    //            result += aij * X[i_col];
    //
    //            if (i_col!=i_row) Y[i_col] += aij * (*x_iterator); 
    //
    //            data_iterator++;
    //        }
    //        *y_iterator += result;
    //
    //        row_iterator++;
    //        y_iterator++;
    //        x_iterator++;
    //        i_row++;
    //    }  
    //}

    //*************************************************************************************************
    //*  Y = A*X with A symmetric matrix - Paralelized OpenMP
    //*************************************************************************************************
    void cMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;
    
        int i_col = 0;
        
		ComplexType aij;
		ComplexType result;
        ComplexType czero(0.0,0.0);

        std::fill(Y.begin(), Y.end(), ComplexType());
		
		omp_set_num_threads(NumThreads);

		int nthreads = omp_get_num_threads();

		#pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
		{
            ComplexVectorType TC(Y);

			#pragma omp for 
			for (int i = 0; i < A.RowsNumber(); i++)
			{
			    data_iterator     = A[i].begin();
				end_data_iterator = A[i].end();

				result = czero;

				while(data_iterator != end_data_iterator)
				{
					i_col = data_iterator->first;
					aij   = data_iterator->second;

					result += aij * X[i_col];
                    
                    if (i_col != i) TC[i_col] += aij * X[i];

					data_iterator++;
				}

				TC[i] += result; 
            }	

            #pragma omp critical  
			for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
	}

	//*************************************************************************************************
	//*  Y = A*X with A Hermitian matrix ( aij = conj(aji) ) - Paralelized OpenMP
	//*************************************************************************************************
	void cMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
	{
		ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
		ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

		int i_col = 0;

		ComplexType aij;
		ComplexType result;
		ComplexType czero(0.0, 0.0);

		std::fill(Y.begin(), Y.end(), ComplexType());

		omp_set_num_threads(NumThreads);

		int nthreads = omp_get_num_threads();

        #pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
		{
			ComplexVectorType TC(Y);

            #pragma omp for 
			for (int i = 0; i < A.RowsNumber(); i++)
			{
				data_iterator     = A[i].begin();
				end_data_iterator = A[i].end();

				result = czero;

				while (data_iterator != end_data_iterator)
				{
					i_col = data_iterator->first;
					aij   = data_iterator->second;

					result += aij * X[i_col];

					if (i_col != i) TC[i_col] += std::conj(aij) * X[i];

					data_iterator++;
				}

				TC[i] += result;
			}

            #pragma omp critical  
			for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
		}
	}

    //*************************************************************************************************
    //*  Y = A * X with A full non-symmetric sparse - Paralelized OpenMP
    //*************************************************************************************************
    void cMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        int i_col = 0;

        ComplexType aij;
        ComplexType result;
        ComplexType czero(0.0, 0.0);

        std::fill(Y.begin(), Y.end(), ComplexType());

        omp_set_num_threads(NumThreads);

        int nthreads = omp_get_num_threads();

        #pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
        {
            ComplexVectorType TC(Y);

            #pragma omp for 
            for (int i = 0; i < A.RowsNumber(); i++)
            {
                data_iterator     = A[i].begin();
                end_data_iterator = A[i].end();

                result = czero;

                while (data_iterator != end_data_iterator)
                {
                    i_col = data_iterator->first;
                    aij   = data_iterator->second;

                    result += aij * X[i_col];

                    data_iterator++;
                }

                TC[i] += result;
            }

            #pragma omp critical  
            for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
    }
    
    //*******************************************************************************************************
    //*  Y = Atrans * X  with A full non-symmetric sparse - Paralelized OpenMP
    //*******************************************************************************************************
    void cTransMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        ComplexType aij;
        int i_col;

        std::fill(Y.begin(), Y.end(), ComplexType());

        omp_set_num_threads(NumThreads);

        int nthreads = omp_get_num_threads();

        #pragma omp parallel  private( data_iterator, end_data_iterator, i_col, aij ) shared (Y, A, X)   
        {
            ComplexVectorType TC(Y);

            #pragma omp for 
            for (int i = 0; i < A.RowsNumber(); i++)
            {
                data_iterator     = A[i].begin();
                end_data_iterator = A[i].end();

                while (data_iterator != end_data_iterator)
                {
                    aij   = data_iterator->second;
                    i_col = data_iterator->first;

                    TC[i_col] += aij * X[i];

                    data_iterator++;
                }
            }

            #pragma omp critical  
            for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
    }

    //*******************************************************************************************************
    //*  Y = Atrans * X  with A hermitic - Paralelized OpenMP
    //*******************************************************************************************************
    void cTransMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        ComplexType aij;
        int i_col;

        std::fill(Y.begin(), Y.end(), ComplexType());

        omp_set_num_threads(NumThreads);

        int nthreads = omp_get_num_threads();

        #pragma omp parallel  private( data_iterator, end_data_iterator, i_col, aij ) shared (Y, A, X)   
        {
            ComplexVectorType TC(Y);

            #pragma omp for 
            for (int i = 0; i < A.RowsNumber(); i++)
            {
                data_iterator     = A[i].begin();
                end_data_iterator = A[i].end();

                while (data_iterator != end_data_iterator)
                {
                    aij   = data_iterator->second;
                    i_col = data_iterator->first;

                    TC[i_col] += aij * X[i];

                    if (i_col != i) TC[i] += std::conj(aij) * X[i_col];

                    data_iterator++;
                }
            }

            #pragma omp critical  
            for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
    }

    //*******************************************************************************************************
    //*  Y = Atrans * X  with A symmetric - Paralelized OpenMP ( same as cMultSymm() )
    //*******************************************************************************************************
    void cTransMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        int i_col = 0;

        ComplexType aij;
        ComplexType result;
        ComplexType czero(0.0,0.0);

        std::fill(Y.begin(), Y.end(), ComplexType());

        omp_set_num_threads(NumThreads);

        int nthreads = omp_get_num_threads();

        #pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
        {
            ComplexVectorType TC(Y);

            #pragma omp for 
            for (int i = 0; i < A.RowsNumber(); i++)
            {
                data_iterator     = A[i].begin();
                end_data_iterator = A[i].end();

                result = czero;

                while(data_iterator != end_data_iterator)
                {
                    i_col = data_iterator->first;
                    aij   = data_iterator->second;

                    result += aij * X[i_col];

                    if (i_col != i) TC[i_col] += aij * X[i];

                    data_iterator++;
                }

                TC[i] += result; 
            }	

            #pragma omp critical  
            for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
    }

    //******************************************************************************************************
    //*   Y = conj(A) * X with A full non-symmetric sparse - Paralelized OpenMP
    //******************************************************************************************************
    void cConjgMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        int i_col = 0;

        ComplexType aij;
        ComplexType result;
        ComplexType czero(0.0,0.0);

        std::fill(Y.begin(), Y.end(), ComplexType());

        omp_set_num_threads(NumThreads);

        int nthreads = omp_get_num_threads();

        #pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
        {
            ComplexVectorType TC(Y);

            #pragma omp for 
            for (int i = 0; i < A.RowsNumber(); i++)
            {
                data_iterator     = A[i].begin();
                end_data_iterator = A[i].end();

                result = czero;

                while(data_iterator != end_data_iterator)
                {
                    i_col = data_iterator->first;
                    aij   = std::conj(data_iterator->second);

                    result += aij * X[i_col];

                    data_iterator++;
                }

                TC[i] += result; 
            }	

            #pragma omp critical  
            for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
    }

	////*************************************************************************************
	////*  Y = conj(A)*X with A symmetric matrix
	////*************************************************************************************
	//void cConjgMultSym( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
	//{
	//    ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
	//    ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();
    //
	//    ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
	//    ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;
    //
	//    ComplexVectorType::iterator y_iterator = Y.begin();
	//    ComplexVectorType::iterator x_iterator = X.begin();
    //
	//    std::fill(Y.begin(), Y.end(), ComplexType());
	//       
	//    ComplexType aij;
    //
	//    int i_row = 0;
	//    int i_col = 0;
    //
	//    while(row_iterator != end_iterator)
	//    {
	//        data_iterator     = row_iterator->begin();
	//        end_data_iterator = row_iterator->end();
    //
	//        ComplexType result(0.00,0.00);
    //
	//        while(data_iterator != end_data_iterator)
	//        {
	//            i_col = data_iterator->first;
	//            aij   = std::conj(data_iterator->second);
    //
	//            result += aij * X[i_col];
    //
	//            if (i_col!=i_row) Y[i_col] += aij * (*x_iterator); 
    //
	//            data_iterator++;
	//        }
	//        *y_iterator += result;
    //
	//        row_iterator++;
	//        y_iterator++;
	//        x_iterator++;
	//        i_row++;
	//    }
	//}

    //********************************************************************************************************
    //*   Y = conj(A)*X with A symmetric matrix - Paralelized OpenMP
    //********************************************************************************************************
    void cConjgMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;
    
        int i_col = 0;
        
		ComplexType aij;
		ComplexType result;
        ComplexType czero(0.0,0.0);

        std::fill(Y.begin(), Y.end(), ComplexType());
		
		omp_set_num_threads(NumThreads);

		int nthreads = omp_get_num_threads();

		#pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
		{
            ComplexVectorType TC(Y);

			#pragma omp for 
			for (int i = 0; i < A.RowsNumber(); i++)
			{
			    data_iterator     = A[i].begin();
				end_data_iterator = A[i].end();

				result = czero;

				while(data_iterator != end_data_iterator)
				{
					i_col = data_iterator->first;
					aij   = std::conj(data_iterator->second);

					result += aij * X[i_col];

                    if (i_col != i) TC[i_col] += aij * X[i];

					data_iterator++;
				}

				TC[i] += result; 
            }	

            #pragma omp critical  
			for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
        }
	}

	//******************************************************************************************************
	//*   Y = conj(A)*X with A Hermitian matrix ( aij = conj(aji) ) - Paralelized OpenMP
	//******************************************************************************************************
	void cConjgMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
	{
		ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
		ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

		int i_col = 0;

		ComplexType aij;
		ComplexType result;
		ComplexType czero(0.0, 0.0);

		std::fill(Y.begin(), Y.end(), ComplexType());

		omp_set_num_threads(NumThreads);

		int nthreads = omp_get_num_threads();

        #pragma omp parallel  private(data_iterator, end_data_iterator, i_col, aij, czero, result) shared (Y, A, X)   
		{
			ComplexVectorType TC(Y);

            #pragma omp for 
			for (int i = 0; i < A.RowsNumber(); i++)
			{
				data_iterator     = A[i].begin();
				end_data_iterator = A[i].end();

				result = czero;

				while (data_iterator != end_data_iterator)
				{
					i_col = data_iterator->first;
					aij   = std::conj(data_iterator->second);

					result += aij * X[i_col];

					if (i_col != i) TC[i_col] += std::conj(aij) * X[i];

					data_iterator++;
				}

				TC[i] += result;
			}

            #pragma omp critical  
			for (int j = 0; j < Y.size(); j++) Y[j] += TC[j];
		}
	}

    //*****************************************************************************
    //*  Y = A*X
    //*****************************************************************************
    void cMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        ComplexVectorType::iterator c_iterator = Y.begin();

        while(row_iterator != end_iterator)
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            ComplexType result(0.00,0.00);

            while(data_iterator != end_data_iterator)
            {
                result += data_iterator->second * X[data_iterator->first];
                data_iterator++;
            }

            *c_iterator = result;

            row_iterator++;
            c_iterator++;
        }
    }

    //**********************************************************************************
    //*  Y = conj(A)*X
    //**********************************************************************************
    void cConjgMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        ComplexVectorType::iterator c_iterator = Y.begin();

        while(row_iterator != end_iterator)
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            ComplexType result(0.00,0.00);

            while(data_iterator != end_data_iterator)
            {
                result += std::conj(data_iterator->second) * X[data_iterator->first]; 
                data_iterator++;
            }

            *c_iterator = result;

            row_iterator++;
            c_iterator++;
        }
    }

    //**************************************************************************************
    //*  Y = At * X
    //**************************************************************************************
    void cTransposeMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        std::fill(Y.begin(), Y.end(), ComplexType());

        ComplexVectorType::iterator c_iterator = X.begin();

        while(row_iterator != end_iterator)
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            while(data_iterator != end_data_iterator)
            {
                Y[data_iterator->first] += data_iterator->second * *c_iterator;
                data_iterator++;
            }
            row_iterator++;
            c_iterator++;
        }
    }

    //*******************************************************************************************
    //*  Y = conj(At) * X
    //*******************************************************************************************
    void cTransposeConjgMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        std::fill(Y.begin(), Y.end(), ComplexType());

        ComplexVectorType::iterator c_iterator = X.begin();

        while(row_iterator != end_iterator)
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            while(data_iterator != end_data_iterator)
            {
                Y[data_iterator->first] += std::conj(data_iterator->second) * *c_iterator;
                data_iterator++;
            }
            row_iterator++;
            c_iterator++;
        }
     }

    //*****************************************
    //* ||X|| = X.real()^2 + X.imag()^2
    //*****************************************
    double cTwoNorm( ComplexVectorType& X )
    {
        ComplexVectorType::iterator x_it  = X.begin();
        ComplexVectorType::iterator x_end = X.end();

        double result = 0.00;

        while (x_it != x_end)
        {
            result += std::norm(*x_it);
            x_it++;
        }
        return result;
    } 

    //**************************************************************
    //*  X*Y = X[0]*Y[0] + X[1]*Y[1] + ... 
    //**************************************************************
    ComplexType cDot( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
         ComplexVectorType::iterator x_it  = X.begin();
         ComplexVectorType::iterator x_end = X.end();
         ComplexVectorType::iterator y_it  = Y.begin();

         ComplexType result(0.00,0.00);

         while (x_it != x_end)
         {
             result += (*x_it) * (*y_it);
             y_it++;
             x_it++;
         }
         return result;
    } 

    //******************************************************************
    //*  X*Y = X[0]*conj(Y[0]) + X[1]*conj(Y[1]) + ... 
    //******************************************************************
    ComplexType cConjgDot( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
         ComplexVectorType::iterator x_it  = X.begin();
         ComplexVectorType::iterator x_end = X.end();
         ComplexVectorType::iterator y_it  = Y.begin();

         ComplexType result(0.00,0.00);

         while (x_it != x_end)
         {
             result += (*x_it) * std::conj(*y_it);
             y_it++;
             x_it++;
         }
         return result;
    } 

	//*****************************************************
    //*  Y =  X + Y 
    //*****************************************************
    void cAdd( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = (*x_iterator) + (*y_iterator);
             y_iterator++;
             x_iterator++;
        }
    } 

    //***************************************************************************
    //*  Z =  X + Y 
    //***************************************************************************
    void cAdd( ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z ) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator z_iterator   = Z.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
            *z_iterator = (*x_iterator) + (*y_iterator);
            x_iterator++;
            y_iterator++;
            z_iterator++;
        }
    } 

	//********************************************************
    //*  Y =  X - Y 
    //********************************************************
    void cDeduct( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = (*x_iterator) - (*y_iterator);
             y_iterator++;
             x_iterator++;
        }
    } 

    //*******************************************************************************
    //*  Z =  X - Y 
    //*******************************************************************************
    void cDeduct( ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z ) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator z_iterator   = Z.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
            *z_iterator = (*x_iterator) - (*y_iterator);
            x_iterator++;
            y_iterator++;
            z_iterator++;
        }
    } 
}