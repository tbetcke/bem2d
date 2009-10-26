#ifndef _INCIDENTFUN_H_
#define _INCIDENTFUN_H_

#include "bem2ddefs.h"
#include "Point.h"
#include "boost/shared_ptr.hpp"
#include <complex>

namespace bem2d {
	
	class Incidentfun {
	public:
		
		virtual ~Incidentfun() {};
		virtual complex operator()(Point p) const=0;
		void setnormal(Point normal);
	protected:
		Point n;
	};
	
	class PlaneWave: public Incidentfun {
	public:
		PlaneWave(Point direction, freqtype kvalue);
		
		inline complex operator()(Point p) const {
			complex i(0,1);
			return std::exp(k*i*(dir.x*p.x+dir.y*p.y));
		}				
	private:
		Point dir;
		freqtype k;
	};
	
	class NormalPlaneWave: public Incidentfun {
	public:
		NormalPlaneWave(Point direction, freqtype kvalue);
		
		inline complex operator()(Point p) const {
			complex i(0,1);
			complex f=std::exp(k*i*(dir.x*p.x+dir.y*p.y));
			return i*k*f*(n.x*p.x+n.y*p.y);
		}
	private:
		Point dir;
		freqtype k;
	};
		
	
	
	typedef boost::shared_ptr<Incidentfun> pIncidentfun;
}
	
#endif // _INCIDENTFUN_H_
