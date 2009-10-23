#include "quadrature.h"
#include <tr1/array>
#include <cstdlib>
#include<boost/shared_ptr.hpp>

namespace bem2d {
	
	
	void gauss(dvector& x, dvector& w, int N) throw (lapack_error) {
		// Points xa and weights w for N point Gaussian quadrature
		// in the interval [0,1]
		
		
        char c = 'V';
        int info;
		
		
        x.clear(); w.clear(); // Make sure that x and w are empty
        x.resize(N); w.resize(N);
		
        if (N == 1) {
            x[0]=0.5;
            w[0]=1.0;
            return;
        }
		
        dvector alpha(N,0);
        dvector beta(N-1,0);
        dvector v(N*N,0);
        dvector work(2*N-2,0);
		
		
        /* Fill the Array */
		
        for (std::size_t i = 0; i < N - 1; i++) {
            beta[i] = .5 / sqrt(1. - 1. / (4. * (i + 1)*(i + 1)));
        }
		
        /* Call Lapack */
		
        dstev_(&c, &N, &alpha[0], &beta[0], &v[0], &N, &work[0], &info);
		
        /* Return results */
		
        for (std::size_t i = 0; i < N; i++) {
            x[i] = (1.0 + alpha[i]) / 2.0;
            w[i] = v[i * N] * v[i * N];
        }
		
        if (info) throw lapack_error();
		
    }
	
	void mappoints(dvector& x, dvector& w, double a, double b){
		
		for (std::size_t i = 0; i < x.size(); i++) {
			x[i] = a + (b - a) * x[i];
			w[i] = (b - a) * w[i];
		}
		
	}
	
	void mappoints2d(dvector& x, dvector& y, dvector& w, double a, double b, double c, double d){
		
		for (std::size_t i=0; i< x.size(); i++){
			x[i]=a+(b-a)*x[i];
			y[i]=c+(d-c)*y[i];
			w[i]=(b-a)*(d-c)*w[i];
		}
	}
	
	complex integrateident(std::pair<pElement,pBasis> testfun, std::pair<pElement,pBasis> basfun, kernel& g, AdaptedGauss1& g2d){
		pElement elem1=testfun.first;
		pElement elem2=basfun.first;
		pBasis fun1=testfun.second;
		pBasis fun2=basfun.second;
		const dvector x=g2d.getx();
		const dvector y=g2d.gety();
		const dvector w=g2d.getw();
		
		complex result=0;
		for (int i=0;i<g2d.size();i++){
			g.setnormal(elem1->normal(x[i]),elem2->normal(y[i]));
			
			complex f1=std::conj((*fun1)(x[i]));
			complex f2=(*fun2)((1-y[i])*x[i]);
			Point xp=elem1->map(x[i]); Point yp=elem2->map((1-y[i])*x[i]);
			complex gvalue=g(xp,yp);
			double s1=length(elem1->deriv(x[i])); double s2=length(elem2->deriv((1-y[i])*x[i]));
			result+= x[i]*gvalue*f1*f2*s1*s2*w[i];
			
			double tau=1-x[i];
			double z=1-x[i]+y[i]*x[i];
			f1=std::conj((*fun1)(tau));			
			f2=(*fun2)(z);
			Point xpp=elem1->map(tau); Point ypp=elem2->map(z);
			gvalue=g(xpp,ypp);
			s1=length(elem1->deriv(tau)); s2=length(elem2->deriv(z));
			result+= tau*gvalue*f1*f2*s1*s2*w[i];
			
			
		}
		return result;
		
		
	}
	
	void evalkernel(std::pair<pElement,pBasis>& basfun, std::vector<Point>& points, 
					cvector& vals, complex alpha, kernel& g, Gauss1D& g1d) throw (size_error) {
		pElement elem=basfun.first;
		pBasis bf=basfun.second;
		const dvector x=g1d.getx();
		const dvector w=g1d.getw();
		if (points.size()!=vals.size()) throw size_error();
		int N=points.size();
		
		for (int j=0;j<N;j++){
			for (int i=0;i<g1d.size();i++){
				g.setnormal(elem->normal(x[i]),elem->normal(x[i]));
				complex basval=(*bf)(x[i]);
				Point xp=elem->map(x[i]);
				double s=length(elem->deriv(x[i]));
				vals[j]+=alpha*g(points[j],xp)*s*basval*w[i];
			}
		}
		
		
	}
	
	
	boost::shared_ptr<std::vector<complex> > discretekernel(const Geometry& geom, QuadOption opts, kernel& g){
		
		boost::shared_ptr<Geometry::flat_basis_map> bfuns=geom.getflatmap();
		boost::shared_ptr<std::vector<complex> > pmatrix(new std::vector<complex >);
		std::size_t N=bfuns->size();
		pmatrix->resize(N*N);
		
		Gauss2D g0(opts.N);
		AdaptedGauss1 g1(opts.N,opts.L,opts.sigma);
		AdaptedGauss2 g2(opts.N,opts.L,opts.sigma);
		AdaptedGauss3 g3(opts.N,opts.L,opts.sigma);
		
		for (int j=0;j<N;j++){
			for (int i=0;i<N;i++){
				std::pair<pElement,pBasis> pi((*bfuns)[i]);
				std::pair<pElement,pBasis> pj((*bfuns)[j]);
				std::size_t indi=pi.first->getIndex();
				std::size_t indj=pj.first->getIndex();
				//std::cout << indi << " " << indj << std::endl;
				
				
				if (indi==indj){
					// Identical elements
					//std::cout << "Integrate Ident" << std::endl;
					(*pmatrix)[i+N*j]=integrateident(pi,pj,g,g1);
				}
				else if (pi.first->getNext()==indj){
					// (1,0) situation
					//std::cout << "Integrate (1,0)" << std::endl;
					(*pmatrix)[i+N*j]=integrate(pi,pj,g,g3);
				}
				else if (pi.first->getPrev()==indj){
					// (0,1) situation
					//std::cout << "Integrate (0,1)" << std::endl;
					(*pmatrix)[i+N*j]=integrate(pi,pj,g,g2);
				}
				else {
					// Elements are not neighbors
					//std::cout << "Integrate Remote" << std::endl;
					(*pmatrix)[i+N*j]=integrate(pi,pj,g,g0);
				}
				//std::cout << (*pmatrix)[i+N*j] << std::endl;
				
			}
		}
		return pmatrix;
		
	}
	
	
	
	
	
}
