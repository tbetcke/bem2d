#ifndef MATHROUTINES_H_
#define MATHROUTINES_H_

#include<vector>
#include "bem2ddefs.h"
#include "exceptions.h"
#include<vector>
#include "Point.h"


namespace bem2d {
	
	complex besselH0(double x);
	complex besselH1(double x);
	
	void abrange(dvector& xx, double a, double b, int n);
	
    // Create a meshgrid, x and y must be of length xpts*ypts
    // ax, bx and ay, by give the grid boundaries and xps, ypts the number of
    // points.
		
	boost::shared_ptr<std::vector<Point > > meshgrid(double ax, double bx,
				  double ay, double by, int xpts, int ypts);
	
	void solve_system(pcvector pmatrix, pcvector prhs) throw (lapack_error);
	
	
}

#endif