#include "bem2d_shape.h"
#include "bem2d_defs.h"

namespace bem2d
{

DiskShapePiecewiseConst::DiskShapePiecewiseConst(int n,double radius): elements_(n), points_(n)
{

  double angle=2*PI/n;
  for (int i=0; i<n; i++) points_[i]=boost::shared_ptr<bem2d::Point>(new bem2d::Point(radius*cos(i*angle),radius*sin(i*angle)));
  for (int i=0; i<n-1; i++)
    {
      elements_[i]=bem2d::pElement(new bem2d::ConstElement(*(points_[i]),*(points_[i+1]),i));
    }
  elements_[n-1]=bem2d::pElement(new bem2d::ConstElement(*(points_[n-1]),*(points_[0]),n-1));
  for (int i=1; i<n-1; i++)
    {
      elements_[i]->set_next(i+1);
      elements_[i]->set_prev(i-1);
    }
  elements_[0]->set_next(1);
  elements_[n-1]->set_next(0);

}

pGeometry DiskShapePiecewiseConst::GetGeometry()
{
  return bem2d::pGeometry(new bem2d::Geometry(elements_));
}

}