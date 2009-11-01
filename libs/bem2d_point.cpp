#include <cmath>
#include "bem2d_point.h"



namespace bem2d {

    Point::Point(double xx, double yy)
    : x(xx),
    y(yy) {
    }

    Point::Point()
    : x(0),
    y(0) {
    }

    Point::Point(const Point& p)
    : x(p.x),
    y(p.y) {
    }


    std::ostream & operator<<(std::ostream &s, const Point& z) {
        return s << "(" << z.x << "," << z.y << ")";
    }

    Point operator+(const Point& p1, const Point& p2){
        return Point(p1.x+p2.x,p1.y+p2.y);
    }

    Point operator-(const Point& p1, const Point& p2){
        return Point(p1.x-p2.x,p1.y-p2.y);
    }

    Point operator*(const Point& p, double t){
        return Point(t*p.x,t*p.y);
    }
    Point operator*(double t,const Point& p){
        return Point(t*p.x,t*p.y);
    }

    double length(const Point& p){
        return sqrt(p.x*p.x+p.y*p.y);
    }

    Point normalize(const Point& p){
        double l=length(p);
        return Point(p.x/l,p.y/l);
    }

}
