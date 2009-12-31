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

template<typename T>
class AnalyticCurve
{
public:
        AnalyticCurve(int n, const T& curve);
        pGeometry GetGeometry();
private:
        std::vector<bem2d::pElement> elements_;
};

template<typename T>
AnalyticCurve<T>::AnalyticCurve(int n, const T& curve): elements_(n)
{

        double h=1.0/n;
        for (int i=0; i<n-1; i++) {
                elements_[i]=pElement(new AnalyticCurveElement<T>(i*h, (i+1)*h,curve,i));
        }
        elements_[n-1]=pElement(new AnalyticCurveElement<T>((n-1)*h, 1,curve,n-1));
        for (int i=1; i<n-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(n-1);
        elements_[n-1]->set_next(0);
        elements_[n-1]->set_prev(n-2);
}

template<typename T>
pGeometry AnalyticCurve<T>::GetGeometry()
{
        return pGeometry(new bem2d::Geometry(elements_));
}




}

#endif // _SHAPE_H_
