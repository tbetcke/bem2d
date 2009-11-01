#ifndef _BASISTYPES_H
#define	_BASISTYPES_H

#include<cstdlib>
#include <complex>
#include "gsl/gsl_sf_legendre.h"
#include "bem2d_defs.h"
#include "bem2d_basis.h"
#include "bem2d_geometry.h"

namespace bem2d{


	
	class ConstBasis: public Basis {
	public:
		ConstBasis();
		inline complex operator()(double t){
			return 1.0;
		}
		
	};
	
	typedef boost::shared_ptr<ConstBasis> pConstBasis;
	
	class PolBasis: public Basis {
	public:
		static void AddBasis(int deg, pGeometry pgeom);			
		PolBasis(int degree);
		inline complex operator()(double t){
			double result=1;
			for (int i=0;i<degree_;i++) result=result*t;
			return result;
		}
	private:
		int degree_;
	};
	
	typedef boost::shared_ptr<PolBasis> pPolBasis;
	
	class WaveBasis: public Basis {
	public:
		static void AddBasis(freqtype k, pGeometry pgeom);
		WaveBasis(double direction, freqtype k);
		inline complex operator()(double t){
			return std::exp(complex(0,1)*k_*direction_*t);
		}
	private:
		double direction_;
		freqtype k_;
	};
	
	typedef boost::shared_ptr<WaveBasis> pWaveBasis;

	class WavePolBasis: public Basis {
	public:
		static void AddBasis(int degree, freqtype k, pGeometry pgeom);
		WavePolBasis(int degree, double direction, freqtype k);
		inline complex operator()(double t){
			double result=1;
			for (int i=0;i<degree_;i++) result=result*t;
			return std::exp(complex(0,1)*k_*direction_*t)*result;
		}
	private:
		double direction_;
		freqtype k_;
		int degree_;
	};
	
	typedef boost::shared_ptr<WavePolBasis> pWavePolBasis;
			
		
}
#endif	/* _BASIS_H */

