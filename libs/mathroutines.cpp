#include "mathroutines.h"
#include "gsl/gsl_sf_bessel.h"
#include "fstream"

namespace bem2d {
	

	complex besselH0(double x) {
        double fr, fi;
		
        // Evaluate Bessel function
		
        fr = gsl_sf_bessel_J0(x);
        fi = gsl_sf_bessel_Y0(x);
        return complex(fr, fi);
		
    }
	
    complex besselH1(double x) {
        double fr, fi;
		
        // Evaluate Bessel function
		
        fr = gsl_sf_bessel_J1(x);
        fi = gsl_sf_bessel_Y1(x);
        return complex(fr, fi);
		
    }
	
	void abrange(dvector& xx, double a, double b, int n){
		xx.resize(n);
		for (int j = 0; j < n; j++) xx[j] = a + (1.0 / (n - 1)) * j * (b - a);
	}

	boost::shared_ptr<std::vector<Point > > meshgrid(double ax, double bx,
												double ay, double by, int xpts, int ypts){
		
		dvector xx;
		dvector yy;
		
		abrange(xx,ax,bx,xpts);
		abrange(yy,ay,by,ypts);
		
		boost::shared_ptr<std::vector<Point > > pvec(new std::vector<Point>);
		pvec->reserve(xpts*ypts);
		
        // Create the grid
		
        for (int i = 0; i < xpts; i++) {
            for (int j = 0; j < ypts; j++) {
				pvec->push_back(Point(xx[i],yy[j]));
            }
        }
		
		return pvec;
    }
	
	void solve_system(pcvector pmatrix, pcvector prhs) throw (lapack_error){
		
		int N=prhs->size();
		int info;
		char trans ='N';
		int nrhs=1;
		std::vector<int> ipiv(N);
		
		
		zgetrf_(&N, &N, &((*pmatrix)[0]), &N, &(ipiv[0]),&info);
		if (info!=0) throw lapack_error();
		
		zgetrs_(&trans,&N,&nrhs,&((*pmatrix)[0]),&N,&(ipiv[0]),&((*prhs)[0]),&N,&info);
		if (info!=0) throw lapack_error();
		
	}
	
	
	
	
	
}
