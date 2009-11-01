#include "bem2d_shape.h"
#include "bem2d_defs.h"

namespace bem2d {

	diskshape_piecewise_const::diskshape_piecewise_const(int n,double radius): elements(n), points(n){
	
		double angle=2*PI/n;
		for (int i=0;i<n;i++) points[i]=boost::shared_ptr<bem2d::Point>(new bem2d::Point(radius*cos(i*angle),radius*sin(i*angle)));
		for (int i=0;i<n-1;i++) {
			elements[i]=bem2d::pElement(new bem2d::ConstElement(*(points[i]),*(points[i+1]),i));
		}
		elements[n-1]=bem2d::pElement(new bem2d::ConstElement(*(points[n-1]),*(points[0]),n-1));
		for (int i=1;i<n-1;i++) {
			elements[i]->set_next(i+1);
			elements[i]->set_prev(i-1);
		}
		elements[0]->set_next(1); elements[n-1]->set_next(0);
		
	}
	
	pGeometry diskshape_piecewise_const::getGeometry(){
		return bem2d::pGeometry(new bem2d::Geometry(elements));
	}
	
}