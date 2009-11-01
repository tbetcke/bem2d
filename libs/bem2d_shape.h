#ifndef _SHAPE_H_
#define _SHAPE_H_

#include<vector>
#include "boost/shared_ptr.hpp"
#include "bem2d_element.h"
#include "bem2d_geometry.h"
#include "bem2d_point.h"
#include "bem2d_defs.h"

namespace bem2d {
	
	class diskshape_piecewise_const {
	public:
		diskshape_piecewise_const(int n,double radius);
		pGeometry getGeometry();
	private:
		std::vector<bem2d::pElement> elements;
		std::vector<boost::shared_ptr<bem2d::Point> > points;

		
	};
	
	template<typename T>
	class analytic_curve {
	public:
		analytic_curve(int n, const T& curve);
		pGeometry getGeometry();
	private:
		std::vector<bem2d::pElement> elements;
	};
	
	template<typename T>
	analytic_curve<T>::analytic_curve(int n, const T& curve): elements(n){
		
		double h=1.0/n;
		for (int i=0;i<n-1;i++){
			elements[i]=pElement(new AnalyticCurveElement<T>(i*h, (i+1)*h,curve,i));
		}
		elements[n-1]=pElement(new AnalyticCurveElement<T>((n-1)*h, 1,curve,n-1));
		for (int i=1;i<n-1;i++) {
			elements[i]->set_next(i+1);
			elements[i]->set_prev(i-1);
		}
		elements[0]->set_next(1); elements[n-1]->set_next(0);
	}
	
	template<typename T>
	pGeometry analytic_curve<T>::getGeometry(){
		return pGeometry(new bem2d::Geometry(elements));
	}
		
		

		
}

#endif // _SHAPE_H_
