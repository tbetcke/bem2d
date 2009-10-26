#ifndef _KERNEL_H
#define _KERNEL_H

#include <boost/utility.hpp>
#include "Point.h"
#include <complex>
#include "bem2ddefs.h"

namespace bem2d {

	class kernel: private boost::noncopyable {
	public:
		virtual ~kernel(){};
		inline void setnormal(Point normal1, Point normal2){
			n1[0]=normal1.x; n1[1]=normal1.y;
			n2[0]=normal2.x; n2[1]=normal2.y;
		}
		inline void setx(Point p){
			x1=p.x;
			x2=p.y;
		}
		
		virtual complex operator()(Point x, Point y) const=0;
		inline complex operator()(Point y) const {
			return (*this)(Point(x1,x2),y);
		}
		
	protected:
		double n1[2];
		double n2[2];
		double x1; double x2;
	};
	
	class singlelayer: public kernel {
	public:
		singlelayer(freqtype kvalue);
		complex operator ()(Point x, Point y) const;
	private:
		freqtype k;
	};

	class doublelayer: public kernel {
	public:
		doublelayer(freqtype kvalue);
		complex operator ()(Point x, Point y) const;
	private:
		freqtype k;
	};
	
	class conjdoublelayer: public kernel {
	public:
		conjdoublelayer(freqtype kvalue);
		complex operator ()(Point x, Point y) const;
	private:
		freqtype k;
	};
	
		
}


#endif // _KERNEL_H