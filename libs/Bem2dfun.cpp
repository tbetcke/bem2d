#include "Bem2dfun.h"

namespace bem2d {

	
	void Bem2dfun::setnormal(Point p){
		n=p;
	}
	
	Bem2dfun::~Bem2dfun(){}
	
	PlaneWave::PlaneWave(Point direction, freqtype kvalue): dir(direction), k(kvalue){}
	
	NormalPlaneWave::NormalPlaneWave(Point direction, freqtype kvalue):
	dir(direction), k(kvalue) {}
	
	Outwave::Outwave(freqtype kvalue): k(kvalue), s(kvalue){}
	
	complex Outwave::operator()(Point p) const{
		Point zero(0,0);
		return s(zero,p);
	}
	
}