
#if !defined(KRATOS_PROPERTIES_DEFAULT_FUNCTION)
#define KRATOS_PROPERTIES_DEFAULT_FUNCTION

#include "../external_libraries/boost/boost/smart_ptr.hpp"

#include "../includes/variables.h"
#include "../includes/point.h"

namespace Kratos
{
    template<class TDataType>
    class PropertiesDefaultFunction
    {
      public:
       
        typedef boost::shared_ptr<PropertiesDefaultFunction<TDataType> > Pointer;
        
        PropertiesDefaultFunction()
        {
        }
        
        virtual ~PropertiesDefaultFunction()
        {
        }
            
        virtual TDataType operator()( const Variable<TDataType>& rV, const Node& rThisNode, double Time )
        {
            return TDataType();
        }
        
        virtual TDataType operator()( const Variable<TDataType>& rV, const Node& rThisNode )
        {
            return TDataType();
        }
        
        virtual TDataType operator()( const Variable<TDataType>& rV )
        {
            return TDataType();
        }
    }; 
}  

#endif 

