#ifndef _BASIS_H
#define	_BASIS_H

#include<cstdlib>
#include<boost/shared_ptr.hpp>
#include "bem2ddefs.h"
#include "gsl/gsl_sf_legendre.h"

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
	
	class LegendrePolBasis: public Basis {
	public:
		LegendrePolBasis(int degree);
		inline complex operator()(double t){
			//return gsl_sf_legendre_Pl(deg,-1.0+2*t);
			double result=1;
			for (int i=0;i<deg;i++) result=result*t;
			return result;
		}
	private:
		int deg;
	};
		
}
#endif	/* _BASIS_H */

