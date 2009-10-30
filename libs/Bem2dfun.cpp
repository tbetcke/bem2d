#include "Bem2dfun.h"

namespace bem2d {

	
	PlaneWave::PlaneWave(Point direction, freqtype kvalue): dir(direction), k(kvalue){}
	PlaneWave::PlaneWave(const PlaneWave& p): k(p.getk()), dir(p.getdir()){}
	
	NormalPlaneWave::NormalPlaneWave(Point direction, freqtype kvalue):
	dir(direction), k(kvalue) {}
	
	NormalPlaneWave::NormalPlaneWave(const NormalPlaneWave& np):
	dir(np.getdir()), k(np.getk()){}

	CombinedPlaneWave::CombinedPlaneWave(Point direction, freqtype kvalue, double etavalue):
	dir(direction), k(kvalue), eta(etavalue) {}
	
	CombinedPlaneWave::CombinedPlaneWave(Point direction, freqtype kvalue):
	dir(direction), k(kvalue), eta(kvalue) {}
	
	
	CombinedPlaneWave::CombinedPlaneWave(const CombinedPlaneWave& np):
	dir(np.getdir()), k(np.getk()), eta(np.eta){}
	
	
	Outwave::Outwave(freqtype kvalue): k(kvalue), s(kvalue){}
	
	Outwave::Outwave(const Outwave& owave): k(owave.getk()), s(owave.getk()){};
	
	complex Outwave::operator()(Point p) const{
		Point zero(0,0);
		return s(zero,p);
	}
	
}