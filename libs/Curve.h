#ifndef _CURVE_H_
#define _CURVE_H_

#include "Point.h"
#include "bem2ddefs.h"
#include<cmath>

namespace bem2d {

	class circle {
	public:
		inline Point map(double t) const {
			return Point(cos(2*PI*t),sin(2*PI*t));
		}
		inline Point deriv(double t) const {
			return 2*PI*Point(-1.0*sin(2*PI*t),cos(2*PI*t));
		}
	};
	
	class trefoil {
	public:
		inline Point map(double t) const {
			double r=1+0.3*cos(3*2*PI*t);
			return Point(r*cos(2*PI*t),r*sin(2*PI*t));
		}
		inline Point deriv(double t) const {
			double r=1+0.3*cos(3*2*PI*t);
			double rp=-0.3*3*2*PI*sin(3*2*PI*t);
			return Point(rp*cos(2*PI*t)-r*2*PI*sin(2*PI*t),rp*sin(2*PI*t)+r*2*PI*cos(2*PI*t));
		}
	};
		
		
}


#endif // _CURVE_H_