#ifndef _SHAPE_H_
#define _SHAPE_H_

#include<vector>
#include "boost/shared_ptr.hpp"
#include "bem2d_element.h"
#include "bem2d_geometry.h"
#include "bem2d_point.h"
#include "bem2d_defs.h"
#include "bem2d_exceptions.h"

namespace bem2d
{

class DiskShapePiecewiseConst
{
public:
        DiskShapePiecewiseConst(int n,double radius);
        pGeometry GetGeometry();
private:
        std::vector<bem2d::pElement> elements_;
        std::vector<boost::shared_ptr<bem2d::Point> > points_;


};

class Polygon
{
public:
        Polygon(const std::vector<Point>& points, int n);
        pGeometry GetGeometry();
private:
        std::vector<bem2d::pElement> elements_;
};


class AnalyticCurve
{
public:
        AnalyticCurve(int n, pCurve curve);
        pGeometry GetGeometry();
	friend int InvAbsDerivative(double t, const double y[], double f[], void* c);
	void ParameterizeArc(int n);
private:
        std::vector<bem2d::pElement> elements_;
	pCurve curve_;
	dvector arclengthparam; // Stores values of t for parameterization with
                                // respect to arc length.
	
};




}

#endif // _SHAPE_H_

