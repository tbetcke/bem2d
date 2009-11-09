#ifndef _POINT_H
#define	_POINT_H

#include<iostream>
#include "bem2d_defs.h"
#include "bem2d_pnpoly.h"

namespace bem2d
{

struct Point
{
  double x; // x coordinate
  double y; // y coordinate

  Point(double xx, double yy);
  Point();
  Point(const Point& p);

  ~Point() {};
};

std::ostream & operator<<(std::ostream &s, const Point& z);
Point operator+(const Point& p1, const Point& p2);
Point operator-(const Point& p1, const Point& p2);
Point operator*(const Point& p, double t);
Point operator*(double t,const Point& p);

double length(const Point& p);
// Return the norm of a point

Point normalize(const Point& p);
// Normalize a point by dividing it through its length

	
	typedef std::vector<Point> PointVector;
	typedef boost::shared_ptr<PointVector> pPointVector;  
}



#endif	/* _POINT_H */

