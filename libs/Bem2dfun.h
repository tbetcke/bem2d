#ifndef _BEM2DFUN_H_
#define _BEM2DFUN_H_

#include "bem2ddefs.h"
#include "Point.h"
#include "boost/shared_ptr.hpp"
#include <complex>
#include "kernel.h"

namespace bem2d {
	
	class Bem2dfun {
	public:
		
		virtual ~Bem2dfun();
		virtual complex operator()(Point p) const=0;
		void setnormal(Point normal);
	protected:
		Point n;
	};
	
	class PlaneWave: public Bem2dfun {
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
	
	class NormalPlaneWave: public Bem2dfun {
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
		
	class Outwave: public Bem2dfun {
	public:
		Outwave(freqtype kvalue);
		complex operator()(Point x) const;
	private:
		freqtype k;
		singlelayer s;
	};
	
	class Idfun: public Bem2dfun {
	public:
		inline complex operator()(Point p) const {
			return 1.0;
		}
	};
		
	typedef boost::shared_ptr<Bem2dfun> pBem2dfun;
}
	
#endif // _INCIDENTFUN_H_
