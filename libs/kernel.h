#ifndef _KERNEL_H
#define _KERNEL_H

#include <boost/utility.hpp>
#include "Point.h"
#include <complex>
#include "bem2ddefs.h"

namespace bem2d {

	class singlelayer {
	public:
		singlelayer(freqtype kvalue);
		singlelayer(const singlelayer& s);
		complex operator ()(Point x, Point y) const;
		inline void setnormal(Point normal1, Point normal2){
			n1=normal1; n2=normal2;		
		}
		inline freqtype getk() const{
			return k;
		}
		
	private:
		freqtype k;
		Point n1; Point n2;
	};

	class doublelayer {
	public:
		doublelayer(freqtype kvalue);
		doublelayer(const doublelayer& d);
		complex operator ()(Point x, Point y) const;
		inline void setnormal(Point normal1, Point normal2){
			n1=normal1; n2=normal2;		
		}
		inline freqtype getk() const{
			return k;
		}
	private:
		freqtype k;
		Point n1; Point n2;
	};
	
	class conjdoublelayer {
	public:
		conjdoublelayer(freqtype kvalue);
		conjdoublelayer(const conjdoublelayer& d);
		complex operator ()(Point x, Point y) const;
		inline void setnormal(Point normal1, Point normal2){
			n1=normal1; n2=normal2;		
		}
		inline freqtype getk() const{
			return k;
		}
	private:
		freqtype k;
		Point n1; Point n2;
	};
	
	class combinedsingleconjdouble {
	public:
		combinedsingleconjdouble(freqtype kvalue, double etavalue);
		combinedsingleconjdouble(freqtype kvalue);
		combinedsingleconjdouble(const combinedsingleconjdouble& scdl);
		complex operator()(Point x, Point y) const;
		inline void setnormal(Point normal1, Point normal2){
			n1=normal1; n2=normal2;
		}
		inline freqtype getk() const{
			return k;
		}
	private:
		freqtype k;
		Point n1; Point n2;
		double eta;
	};
	
	
}


#endif // _KERNEL_H