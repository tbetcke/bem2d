#ifndef _OUTPUTHANDLER_H_
#define _OUTPUTHANDLER_H_

#include "bem2ddefs.h"
#include<vector>
#include "Point.h"
#include <boost/shared_ptr.hpp>

namespace bem2d {
	
	class Outputhandler {
		
		virtual void operator()(std::vector<Point>, cvector vals)=0;
	};

	typedef boost::shared_ptr<Outputhandler> pOutputhandler;
	
}

#endif // _OUTPUTHANDLER_H_