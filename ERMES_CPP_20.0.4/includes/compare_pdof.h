
#if !defined(DEFINE_KRATOS_COMP_PDOF)
#define DEFINE_KRATOS_COMP_PDOF

#include "../includes/dof.h"

namespace Kratos
{
    struct ComparePDof
    {
        bool operator()( Dof::Pointer a1,  Dof::Pointer a2 )
        {
            if ( a1->NodeIndex() ==  a2->NodeIndex() ) 
            {
                return ( a1->GetVariableKey() < a2->GetVariableKey() ); 
            }
            else 
            {
                return ( a1->NodeIndex() <  a2->NodeIndex() );
            }
        }
    };
}  

#endif 

