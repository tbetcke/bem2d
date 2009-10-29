#ifndef _BASIS_H_
#define _BASIS_H_

#include "bem2ddefs.h"

namespace bem2d {

	class Basis {
	public:
		virtual complex operator()(double t)=0;
		virtual ~Basis();
	};
		
	typedef boost::shared_ptr<Basis> pBasis;

}


#endif // _BASIS_H_