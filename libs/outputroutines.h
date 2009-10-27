#ifndef OUTPUTROUTINES_H_
#define OUTPUTROUTINES_H_

#include<vector>
#include "bem2ddefs.h"
#include "Point.h"


namespace bem2d {

	
	void abrange(dvector& xx, double a, double b, int n);
		
	boost::shared_ptr<std::vector<Point > > meshgrid(double ax, double bx,
													 double ay, double by, int xpts, int ypts);	
	
	void gplotout(std::string name, const std::vector<Point>& points, const dvector& z,
				  int xpts, int ypts);
	
}

#endif

