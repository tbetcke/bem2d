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
			n1=normal1; n2=normal2;		
		}
							
		virtual complex operator()(Point x, Point y) const=0;
		
	protected:
		Point n1;
		Point n2;
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