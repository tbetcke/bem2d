#ifndef _KERNEL_H
#define _KERNEL_H

#include <boost/utility.hpp>
#include "Point.h"
#include <complex>
#include "bem2ddefs.h"

namespace bem2d {

	class kernel: private boost::noncopyable {
	public:
		inline void setnormal(Point normal1, Point normal2){
			n1[0]=normal1.x; n1[1]=normal1.y;
			n2[0]=normal2.x; n2[1]=normal2.y;
		}
		virtual complex operator()(Point x, Point y) const=0;
		
	protected:
		double n1[2];
		double n2[2];
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