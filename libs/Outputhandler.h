#ifndef _OUTPUTHANDLER_H_
#define _OUTPUTHANDLER_H_

#include "bem2ddefs.h"
#include<vector>
#include "Point.h"
#include <boost/shared_ptr.hpp>
#include <string>

namespace bem2d {
	
	class Outputhandler {
	public:
		virtual void operator()(const std::vector<Point>& points, const cvector& vals)=0;
	};

	typedef boost::shared_ptr<Outputhandler> pOutputhandler;
	
	class Gplotoutput: public Outputhandler {
	public:
		Gplotoutput(int xptsvalue, int yptsvalue, std::string name);
		void operator()(const std::vector<Point>& points, const cvector& vals);
	private:
		std::string name;
		int xpts;
		int ypts;
	};
	
}

#endif // _OUTPUTHANDLER_H_