#include "Basistypes.h"

namespace bem2d {
	
	ConstBasis::ConstBasis(){
	};

	PolBasis::PolBasis(int degree): deg(degree){};
	
	void PolBasis::addBasis(int deg, pGeometry pgeom){
		
		for (int i=0;i<=deg;i++){
			pgeom->addBasis(pBasis(new PolBasis(i)));
		}
	}
	
	WaveBasis::WaveBasis(double direction, freqtype kvalue): dir(direction), k(kvalue){}
	
	void WaveBasis::addBasis(freqtype k, pGeometry pgeom){
		pgeom->addBasis(pBasis(new WaveBasis(1.0, k)));
		pgeom->addBasis(pBasis(new WaveBasis(-1.0, k)));
	}

	WavePolBasis::WavePolBasis(int degree, double direction, freqtype kvalue): deg(degree), dir(direction), k(kvalue){}
	
	void WavePolBasis::addBasis(int degree, freqtype k, pGeometry pgeom){
		for (int i=0;i<=degree;i++){
			pgeom->addBasis(pBasis(new WavePolBasis(i,1.0, k)));
			pgeom->addBasis(pBasis(new WavePolBasis(i,-1.0, k)));
		}
	}
	
	

}
