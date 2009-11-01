#ifndef OUTPUTROUTINES_H_
#define OUTPUTROUTINES_H_

#include<vector>
#include "bem2d_defs.h"
#include "bem2d_point.h"


namespace bem2d
{


void Abrange(dvector& xx, double a, double b, int n);

boost::shared_ptr<std::vector<Point > > MeshGrid(double ax, double bx,
    double ay, double by, int xpts, int ypts);

void GplotOut(std::string name, const std::vector<Point>& points, const dvector& z,
              int xpts, int ypts);

}

#endif

