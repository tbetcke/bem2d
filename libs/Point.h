/* 
 * File:   Point.h
 * Author: tbetcke
 *
 * Created on July 30, 2009, 7:01 PM
 */
#ifndef _POINT_H
#define	_POINT_H

#include<iostream>

namespace bem2d {

    struct Point {
        const double x; // x coordinate
        const double y; // y coordinate

        Point(double xx, double yy);
        Point();
        Point(const Point& p);

        ~Point(){};
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

}



#endif	/* _POINT_H */

