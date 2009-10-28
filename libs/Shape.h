#ifndef _SHAPE_H_
#define _SHAPE_H_

#include<vector>
#include "Element.h"
#include "geometry.h"
#include <boost/shared_ptr.hpp>

namespace bem2d {
	
	class diskshape_piecewise_const {
	public:
		diskshape_piecewise_const(int n,double radius);
		pGeometry getGeometry();
	private:
		std::vector<bem2d::pElement> elements;
		std::vector<boost::shared_ptr<bem2d::Point> > points;

		
	};
	
}

#endif // _SHAPE_H_
