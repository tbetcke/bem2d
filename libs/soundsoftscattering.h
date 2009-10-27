#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include "bem2ddefs.h"
#include "Geometry.h"
#include "Bem2dfun.h"
#include <vector>
#include "Outputhandler.h"
#include "quadrature.h"
#include "mathroutines.h"

namespace bem2d {
	
	class Soundsoftscattering {
	public:
		Soundsoftscattering(pGeometry geom, freqtype kvalue);
		virtual ~Soundsoftscattering();
	
		void setincoming(pBem2dfun fun);
		void setincoming(pBem2dfun fun, pBem2dfun normalderiv);
		
		
		void setquadoption(int L, int N, double sigma);
		virtual void discretize();
		
		void solve();
		
		pcvector evalincident(const std::vector<Point>& points);
		pcvector evalsol(const std::vector<Point>& points);
		
		void setOutput(pOutputhandler output);
		void plotIncident(const std::vector<Point>& points);
		void plotScattered(const std::vector<Point>& points);
		void plotFull(const std::vector<Point>& points);
		void plotAll(const std::vector<Point>& points);
		
	private:
		pGeometry pgeom;
		freqtype k;
		
		pBem2dfun incoming;
		pBem2dfun incomingnormal;
		
		pOutputhandler pout;
		
		QuadOption quadopts;
		matrix A;
		pcvector prhs;
		pcvector psol;
		
		
	};
	
}


#endif // _SOUNDSOFTSCATTERING_H_
