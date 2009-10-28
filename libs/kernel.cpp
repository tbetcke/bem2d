#include "kernel.h"
#include "mathroutines.h"

namespace bem2d {

	singlelayer::singlelayer(freqtype kvalue){
		k=kvalue;
	}
	
	singlelayer::singlelayer(const singlelayer& s){
		k=s.getk();
	}
	
	complex singlelayer::operator()(Point x, Point y) const{
		
		double d=length(x-y);
		complex i(0,1);
		return i/4.0*(besselH0(k*d));
		
		
	}

	doublelayer::doublelayer(freqtype kvalue){
		k=kvalue;
	}

	doublelayer::doublelayer(const doublelayer& d){
		k=d.getk();
	}
	
	
	complex doublelayer::operator()(Point x, Point y) const {
		Point w=y-x;
		double d=length(w);		
		complex i(0,1);
		
		return -i*k/4.0*(besselH1(k*d)/d*(n2.x*w.x+n2.y*w.y));
	}

	
	
	conjdoublelayer::conjdoublelayer(freqtype kvalue){
		k=kvalue;
	}

	conjdoublelayer::conjdoublelayer(const conjdoublelayer& cd){
		k=cd.getk();
	}
	
	
	complex conjdoublelayer::operator()(Point x, Point y) const {
		Point w=x-y;
		double d=length(w);		
		complex i(0,1);
		
		return -i*k/4.0*(besselH1(k*d)/d*(n1.x*w.x+n1.y*w.y));
	}
}