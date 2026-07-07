
#include <iostream>
#include "complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Initializing private properties
    //*********************************************************************    
    ComplexSolver::ComplexSolver( unsigned int NumThreads, unsigned int MaxIterations, unsigned int StepIterations, double Tolerance, String& Precond )
    {
        mNumThreads      = NumThreads;
		mMaxIterations   = MaxIterations;
		mStepIterations  = StepIterations;
        mTolerance       = Tolerance; 
        mResidualNorm    = 0.0;
		mResEstimation   = 0.0;
        msPreconditioner = Precond;

		mWriteSolutionEveryStep = false;
    }

    //*********************************************************************
    //* - Show info
    //*********************************************************************    
    void ComplexSolver::PrintInfo()
    {
        if (mBNorm == 0.00)
        {
            std::cout << "!!!!!!!!!!!!!!! ||b|| = 0 !!!!!!!!!!!!!!!" << std::endl;
        }
        else if (mResEstimation > mTolerance)
        {
            std::cout << "!!!!!!!!!!! NO CONVERGENCE !!!!!!!!!!!" << std::endl;
			std::cout << std::endl;
            std::cout << "||b||             = " << mBNorm                 << std::endl;
            std::cout << "||Ax-b||          = " << mResidualNorm          << std::endl;
            std::cout << "||Ax-b||/||b||    = " << mResidualNorm / mBNorm << std::endl;
			std::cout << "Solver estimation = " << mResEstimation         << std::endl;
        }
        else
        {
            std::cout << "--------> CONVERGENCE ACHIEVED <-------- " << std::endl;  
            std::cout << std::endl;
            std::cout << "||b||             = " << mBNorm                 << std::endl;
            std::cout << "||Ax-b||          = " << mResidualNorm          << std::endl;
            std::cout << "||Ax-b||/||b||    = " << mResidualNorm / mBNorm << std::endl;
			std::cout << "Solver estimation = " << mResEstimation         << std::endl;
        }

        std::cout << std::endl;

        std::cout << "Max. iterations: " << mMaxIterations << std::endl;
        std::cout << "Tolerance      : " << mTolerance     << std::endl;
        std::cout << "Problem size   : " << mProblemSize   << std::endl;
    
        std::cout << std::endl;

        std::cout << "Iterations     : " << mIterationNumber                      << std::endl;
        std::cout << "Time spent     : " << difftime(mcEnd,mcStart) << " seconds" << std::endl;

        std::cout << std::endl; 
    }

    //***********************************************************
    //* - Initialization of iterative solver
    //***********************************************************
    bool ComplexSolver::IniSolver( char* SolverName, ComplexMatrixType& A, ComplexVectorType& b )
    {
        std::cout << "Iterative solver: " << SolverName       << std::endl;
        std::cout << "Preconditioner  : " << msPreconditioner << std::endl;
        std::cout << std::endl;

        mProblemSize     = b.size();
        mIterationNumber = 0;
        mInfoCount       = 0;
 
        time(&mcStart);

        mBNorm = sqrt(cTwoNorm(b));

        if (mBNorm == 0.00) 
        {
            time(&mcEnd);
            return false;
        }

        return true;
    }

	//***********************************************************
    //* - Write the x_vector every step iteration
    //***********************************************************
    void ComplexSolver::SetWriteStepResult( bool ToWriteStepResult )
    {
		mWriteSolutionEveryStep = ToWriteStepResult;
    }

    //***********************************************************
    //* - Finishing iterative solver
    //***********************************************************
    void ComplexSolver::FinishSolver()
    {
        time(&mcEnd);

        std::cout << std::endl;

        PrintInfo();
    }

}