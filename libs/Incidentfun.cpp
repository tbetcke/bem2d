#include "Incidentfun.h"

namespace bem2d {

	PlaneWave::PlaneWave(Point direction, freqtype kvalue): dir(direction), k(value){}
	
	NormalPlaneWave::NormalPlaneWave(Point direction, freqtype kvalue):
	dir(direction), k(kvalue) {}
	

}