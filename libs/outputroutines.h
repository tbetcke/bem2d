#ifndef OUTPUTROUTINES_H_
#define OUTPUTROUTINES_H_

#include<vector>
#include "bem2ddefs.h"
#include "Point.h"


namespace bem2d {
	
	void gplotout(std::string name, std::vector<Point>& points, dvector& z,
				  int xpts, int ypts);
	
}

#endif

