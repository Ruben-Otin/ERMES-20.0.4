
#if !defined(COMPLEX_MATRIX_OPERATIONS_H_INCLUDED)
#define COMPLEX_MATRIX_OPERATIONS_H_INCLUDED

#include <complex>
#include "../includes/vector.h"
#include "../includes/csr_matrix.h"

#include <math.h>
#include <numeric>

namespace Kratos
{
    typedef std::complex<double> ComplexType;

    typedef CSRMatrix<ComplexType> ComplexMatrixType;

    typedef std::vector<ComplexType> ComplexVectorType;

    void cMultSymm     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
	void cConjgMultSymm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);

	void cMultHerm     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
	void cConjgMultHerm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);

    void cMultFull     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cConjgMultFull(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultFull(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultHerm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultSymm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
  
    void cMult         (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);
    void cConjgMult    (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);

    void cTransposeMult     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);
    void cTransposeConjgMult(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);
    
    double cTwoNorm(ComplexVectorType& X);

    ComplexType cDot     (ComplexVectorType& X, ComplexVectorType& Y);
    ComplexType cConjgDot(ComplexVectorType& X, ComplexVectorType& Y); 

	void cAdd(ComplexVectorType& X, ComplexVectorType& Y); 
    void cAdd(ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z);

	void cDeduct(ComplexVectorType& X, ComplexVectorType& Y); 
    void cDeduct(ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z); 

    //*****************************************
    //* rY = (c * rY) 
    //*****************************************
    template<class TDataType>
    void cScale(TDataType c, ComplexVectorType& Y) 
    {
	    ComplexVectorType::iterator y_it  = Y.begin();
        ComplexVectorType::iterator y_end = Y.end();
        
        while (y_it != y_end)
	    {
	          *y_it = c * (*y_it);
              y_it++;
	    }
    } 

    //*****************************************
    //*  Y = (a * X)
    //*****************************************
    template<class TDataType>
    void cScale(TDataType a, 
                ComplexVectorType& X, 
                ComplexVectorType& Y) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = a * (*x_iterator);
             y_iterator++;
             x_iterator++;
        }
    } 

    //*****************************************
    //*  Z = (a * X) + (b * Y)
    //*****************************************
    template<class TDataType, class TDataType2>
    void cScaleAndAdd(TDataType  a, 
                      ComplexVectorType& X, 
                      TDataType2 b, 
                      ComplexVectorType& Y, 
                      ComplexVectorType& Z)  
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator z_iterator   = Z.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
            *z_iterator++ = (a * *x_iterator++) + (b * *y_iterator++);
        }
    } 

    //***********************************************
    //*  Y = (a * X) + (b * Y) 
    //***********************************************
    template<class TDataType, class TDataType2>
    void cScaleAndAdd(TDataType  a, 
                      ComplexVectorType& X, 
                      TDataType2 b, 
                      ComplexVectorType& Y) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = (a * (*x_iterator)) + (b * (*y_iterator));
             y_iterator++;
             x_iterator++;
        }
    } 

    //***********************************************
    //*  Y = (a * X) +  Y 
    //***********************************************
    template<class TDataType>
    void cScaleFirstAndAdd(TDataType a, 
                           ComplexVectorType& X, 
                           ComplexVectorType& Y) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = (a * (*x_iterator)) + (*y_iterator);
             y_iterator++;
             x_iterator++;
        }
    } 

	//*****************************************
    //*  Z = (a * X) +  Y
    //*****************************************
    template<class TDataType>
    void cScaleFirstAndAdd(TDataType a, 
                           ComplexVectorType& X, 
                           ComplexVectorType& Y, 
                           ComplexVectorType& Z)  
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator z_iterator   = Z.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
            *z_iterator++ = (a * *x_iterator++) + (*y_iterator++);
        }
    } 
 
    //***********************************************
    //*  Y =  X + (b * Y) 
    //***********************************************
    template<class TDataType>
    void cScaleSecondAndAdd(ComplexVectorType& X, 
                            TDataType b, 
                            ComplexVectorType& Y) 
    {
        ComplexVectorType::iterator x_iterator   = X.begin();
        ComplexVectorType::iterator y_iterator   = Y.begin();
        ComplexVectorType::iterator end_iterator = X.end();

        while(x_iterator != end_iterator)
        {
             *y_iterator = (*x_iterator) + (b * (*y_iterator));
             y_iterator++;
             x_iterator++;
        }
    } 
}

#endif