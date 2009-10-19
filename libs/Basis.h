#ifndef _BASIS_H
#define	_BASIS_H

#include<cstdlib>
#include<boost/shared_ptr.hpp>
#include "bem2ddefs.h"

namespace bem2d{
	
	class Basis {
	public:
		virtual complex operator()(double t)=0;
		virtual ~Basis();
	};
		
	typedef boost::shared_ptr<Basis> pBasis;
	
	class ConstBasis: public Basis {
	public:
		ConstBasis();
		inline complex operator()(double t){
			return 1.0;
		}
		
	};
}
#endif	/* _BASIS_H */

