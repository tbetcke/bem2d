#ifndef _BASISTYPES_H
#define	_BASISTYPES_H

#include<cstdlib>
#include<boost/shared_ptr.hpp>
#include "bem2ddefs.h"
#include "gsl/gsl_sf_legendre.h"
#include "Basis.h"
#include "geometry.h"
#include <complex>

namespace bem2d{


	
	class ConstBasis: public Basis {
	public:
		ConstBasis();
		inline complex operator()(double t){
			return 1.0;
		}
		
	};
	
	class PolBasis: public Basis {
	public:
		static void addBasis(int deg, pGeometry pgeom);			
		PolBasis(int degree);
		inline complex operator()(double t){
			double result=1;
			for (int i=0;i<deg;i++) result=result*t;
			return result;
		}
	private:
		int deg;
	};
	
	class WaveBasis: public Basis {
	public:
		static void addBasis(freqtype k, pGeometry pgeom);
		WaveBasis(double direction, freqtype kvalue);
		inline complex operator()(double t){
			return std::exp(complex(0,1)*k*dir*t);
		}
	private:
		double dir;
		freqtype k;
	};

	class WavePolBasis: public Basis {
	public:
		static void addBasis(int degree, freqtype k, pGeometry pgeom);
		WavePolBasis(int degree, double direction, freqtype kvalue);
		inline complex operator()(double t){
			double result=1;
			for (int i=0;i<deg;i++) result=result*t;
			return std::exp(complex(0,1)*k*dir*t)*result;
		}
	private:
		double dir;
		freqtype k;
		int deg;
	};
	
			
		
}
#endif	/* _BASIS_H */

