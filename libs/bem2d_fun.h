#ifndef _BEM2DFUN_H_
#define _BEM2DFUN_H_

#include <complex>
#include "boost/shared_ptr.hpp"
#include "bem2d_defs.h"
#include "bem2d_point.h"
#include "bem2d_kernel.h"

namespace bem2d {
	
	
	class PlaneWave {
	public:
		PlaneWave(){};
		PlaneWave(Point direction, freqtype kvalue);
		PlaneWave(const PlaneWave& p);
		
		inline complex operator()(Point p) const {
			complex i(0,1);
			return std::exp(k*i*(dir.x*p.x+dir.y*p.y));
		}				
		void setnormal(Point normal){};
		inline freqtype getk() const{
			return k;
		}
		inline Point getdir() const{
			return dir;
		}
	private:
		Point dir;
		freqtype k;
	};
	
	class NormalPlaneWave {
	public:
		NormalPlaneWave(){};
		NormalPlaneWave(Point direction, freqtype kvalue);
		NormalPlaneWave(const NormalPlaneWave& np);
		
		inline complex operator()(Point p) const {
			complex i(0,1);
			complex f=std::exp(k*i*(dir.x*p.x+dir.y*p.y));
			return i*k*f*(n.x*dir.x+n.y*dir.y);
		}
		inline void setnormal(Point normal){
			n=normal;
		}
		inline freqtype getk() const{
			return k;
		}
		inline Point getdir() const{
			return dir;
		}
	private:
		Point dir;
		freqtype k;
		Point n;
	};
	
	class CombinedPlaneWave {
	public:
		CombinedPlaneWave(Point direction, freqtype kvalue, double etavalue);
		CombinedPlaneWave(Point direction, freqtype kvalue);
		CombinedPlaneWave(const CombinedPlaneWave& np);
		
		inline complex operator()(Point p) const {
			complex i(0,1);
			complex f=std::exp(k*i*(dir.x*p.x+dir.y*p.y));
			return i*k*f*(n.x*dir.x+n.y*dir.y)+i*eta*f;
		}
		inline void setnormal(Point normal){
			n=normal;
		}
		inline freqtype getk() const{
			return k;
		}
		inline Point getdir() const{
			return dir;
		}
		inline double geteta() const{
			return eta;
		}
	private:
		Point dir;
		freqtype k;
		Point n;
		double eta;
	};
	
	
	
		
	class Outwave {
	public:
		Outwave(freqtype kvalue);
		Outwave(const Outwave& owave);
		complex operator()(Point x) const;
		void setnormal(Point normal){};
		inline freqtype getk() const{
			return k;
		}
	private:
		freqtype k;
		SingleLayer s;
	};
	
	class Idfun {
	public:
		inline complex operator()(Point p) const {
			return 1.0;
		}
		void setnormal(Point normal){};
	};
		
}
	
#endif // _INCIDENTFUN_H_
