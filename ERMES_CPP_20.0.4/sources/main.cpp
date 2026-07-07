
#include <iostream>
#include <complex>

#include "../includes/kernel.h"
#include "../includes/variables.h"

#include "../external_libraries/boost/boost/progress.hpp"

int main( int argc, char *argv[ ], char *envp[ ] )
{
    if( argc < 2 ) 
    {
        throw Kratos::Exception( "main()", std::invalid_argument( "data file needed" ) );
    }

    Kratos::Kernel* p_kernel;

    p_kernel = new Kratos::Kernel;	      

    std::cout << "ERMES analysis started..." << std::endl; 
    std::cout << std::endl;

    p_kernel->Execute( argv[ 1 ] );

    std::cout << std::endl;
	std::cout << "ERMES analysis finished." << std::endl;
    std::cout << std::endl;

    delete p_kernel;

	return 0;
}
