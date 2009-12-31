#include "bem2d_shape.h"
#include "bem2d_defs.h"

namespace bem2d
{

DiskShapePiecewiseConst::DiskShapePiecewiseConst(int n,double radius): elements_(n), points_(n)
{

        double angle=2*PI/n;
        for (int i=0; i<n; i++) points_[i]=boost::shared_ptr<bem2d::Point>(new bem2d::Point(radius*cos(i*angle),radius*sin(i*angle)));
        for (int i=0; i<n-1; i++) {
                elements_[i]=bem2d::pElement(new bem2d::ConstElement(*(points_[i]),*(points_[i+1]),i));
        }
        elements_[n-1]=bem2d::pElement(new bem2d::ConstElement(*(points_[n-1]),*(points_[0]),n-1));
        for (int i=1; i<n-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(n-1);
        elements_[n-1]->set_next(0);
        elements_[n-1]->set_prev(n-2);

}

pGeometry DiskShapePiecewiseConst::GetGeometry()
{
        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}

Polygon::Polygon(const std::vector<Point>& points, int n)
{

        int esize=n*points.size();
        elements_.resize(esize);
        std::vector<Point> p(points);
        p.push_back(p[0]); // Add last element again - makes next for-loop easier
        for (int i=0; i<points.size(); i++) {
                Point direction=p[i+1]-p[i];
                for (int j=0; j<n; j++) {
                        Point p1=p[i]+(1.0*j)/n*direction;
                        Point p2=p[i]+(1.0*(j+1))/n*direction;
                        elements_[n*i+j]=bem2d::pElement(new bem2d::ConstElement(p1,p2,n*i+j));
                }
        }
        for (int i=1; i<elements_.size()-1; i++) {
                elements_[i]->set_next(i+1);
                elements_[i]->set_prev(i-1);
        }
        elements_[0]->set_next(1);
        elements_[0]->set_prev(esize-1);
        elements_[n-1]->set_next(0);
        elements_[esize-1]->set_prev(esize-2);

}


pGeometry Polygon::GetGeometry()
{
        return bem2d::pGeometry(new bem2d::Geometry(elements_));
}


}