#ifndef _QUADRATURE_H
#define	_QUADRATURE_H

#include "bem2ddefs.h"
#include "exceptions.h"
#include "kernel.h"
#include <map>
#include "Element.h"
#include  "Basis.h"
#include "QuadPoints.h"
#include <complex>
#include "Point.h"
#include "geometry.h"
#include "boost/shared_ptr.hpp"

namespace bem2d {

	
	struct QuadOption {
		int L;
		int N;
		double sigma;
	};
	
	
    void gauss(dvector& x, dvector& w, int N) throw (lapack_error);
    void mappoints(dvector& x, dvector& w, double a, double b);
    void mappoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d);

	complex integrateident(std::pair<pElement,pBasis> testfun, std::pair<pElement,pBasis> basfun, kernel& g, AdaptedGauss1& g2d);

	
	template<typename T>
	complex integrate(std::pair<pElement,pBasis>& testfun, std::pair<pElement,pBasis>& basfun, kernel& g, T& g2d){
		
		pElement elem1=testfun.first;
		pElement elem2=basfun.first;
		pBasis fun1=testfun.second;
		pBasis fun2=basfun.second;
		const dvector x=g2d.getx();
		const dvector y=g2d.gety();
		const dvector w=g2d.getw();
		std::size_t t1=elem1->getIndex(); std::size_t t2=elem2->getIndex();
		
		complex result=0;
		for (int i=0;i<g2d.size();i++){
			g.setnormal(elem1->normal(x[i]),elem2->normal(y[i]));
			complex f1=std::conj((*fun1)(x[i]));
			complex f2=(*fun2)(y[i]);
			Point xp=elem1->map(x[i]); Point yp=elem2->map(y[i]);
			complex gvalue=g(xp,yp);
			double s1=length(elem1->deriv(x[i])); double s2=length(elem2->deriv(y[i]));
			result+= gvalue*f1*f2*s1*s2*w[i];
		}
		return result;
				
	}

	void evalkernel(std::pair<pElement,pBasis>& basfun, std::vector<Point>& points, 
					cvector& vals, complex alpha, kernel& g, Gauss1D& g1d) throw (size_error);
	
	template<typename T>
	complex integrate(std::pair<pElement,pBasis> testfun, T& fun, Gauss1D& g1d){
		pElement elem=testfun.first;
		pBasis tf=testfun.second;
		const dvector x=g1d.getx();
		const dvector w=g1d.getw();
		
		complex result=0;
		for (int i=0;i<g1d.size();i++){
			fun.setnormal(elem->normal(x[i]));
			complex tfval=std::conj((*tf)(x[i]));
			Point xp=elem->map(x[i]);
			double s=length(elem->deriv(x[i]));
			result+=fun(xp)*s*tfval*w[i];
		}
		return result;
	}
	
	template<typename T>
	boost::shared_ptr<std::vector<complex> > discreterhs(const Geometry& geom, QuadOption ops, T& fun){
		
		boost::shared_ptr<std::vector<complex> > prhs(new std::vector<complex >);
		int N=geom.getsize();
		prhs->resize(N);
		boost::shared_ptr<Geometry::flat_basis_map> bfuns=geom.getflatmap();
		Gauss1D g1d(ops.N);
	
		for (int i=0;i<N;i++) (*prhs)[i]=integrate((*bfuns)[i],fun, g1d);
		return prhs;
	}
							
	boost::shared_ptr<std::vector<complex> > discretekernel(const Geometry& geom,QuadOption ops, kernel& gl);

	
}
#endif	/* _QUADRATURE_H */

