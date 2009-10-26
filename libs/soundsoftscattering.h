#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include "bem2ddefs.h"
#include "Geometry.h"
#include "Incidentfun.h"
#include <vector>
#include "Outputhandler.h"
#include "quadrature.h"

namespace bem2d {
	
	class Soundsoftscattering {
	public:
		Soundsoftscattering(pGeometry geom, freqtype kvalue);
		~Soundsoftscattering();
	
		void setincoming(pIncidentfun fun);
		void setincoming(pIncidentfun fun, pIncidentfun normalderiv);
		
		
		void setquadoption(int L, int N, double sigma);
		virtual void discretize();
		
		void solve();
		
		void setOutput(pOutputhandler output);
		void plotIncident(std::vector<Point> points);
		void plotScattered(std::vector<Point> points);
		void plotFull(std::vector<Point> points);
		
	private:
		pGeometry pgeom;
		freqtype k;
		
		pIncidentfun incoming;
		pIncidentfun incomingnormal;
		
		pOutputhandler pout;
		
		QuadOption quadopts;
		matrix A;
		pcvector prhs;
		
		
	};
	
}


#endif // _SOUNDSOFTSCATTERING_H_
