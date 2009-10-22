#include "kernel.h"
#include "gsl/gsl_sf_bessel.h"

namespace bem2d {

	singlelayer::singlelayer(freqtype kvalue){
		k=kvalue;
	}
	
	complex singlelayer::operator()(Point x, Point y) const{
		
		double d=length(x-y);
		double bj=gsl_sf_bessel_J0(k*d);
		double by=gsl_sf_bessel_Y0(k*d);

		return complex(0,1)/4.0*(bj+complex(0,1)*by);
		
		
	}

	doublelayer::doublelayer(freqtype kvalue){
		k=kvalue;
	}
	
	complex doublelayer::operator()(Point x, Point y) const {
		Point w=y-x;
		double d=length(w);		
		double bj=gsl_sf_bessel_J1(k*d);
		double by=gsl_sf_bessel_Y1(k*d);
		complex i(0,1);
		
		return -i*k/4.0*((bj+i*by)/d*(n2[0]*w.x+n2[1]*w.y));
	}

	
	
	conjdoublelayer::conjdoublelayer(freqtype kvalue){
		k=kvalue;
	}
	
	complex conjdoublelayer::operator()(Point x, Point y) const {
		Point w=x-y;
		double d=length(w);		
		double bj=gsl_sf_bessel_J1(k*d);
		double by=gsl_sf_bessel_Y1(k*d);
		complex i(0,1);
		
		return -i*k/4.0*((bj+i*by)/d*(n1[0]*w.x+n1[1]*w.y));
	}
}