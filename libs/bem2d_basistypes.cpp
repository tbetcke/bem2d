#include "bem2d_basistypes.h"

namespace bem2d {
	
	ConstBasis::ConstBasis(){
	};

	PolBasis::PolBasis(int degree): degree_(degree){};
	
	void PolBasis::addBasis(int degree, pGeometry pgeom){
		
		for (int i=0;i<=degree;i++){
			pgeom->addBasis(pBasis(new PolBasis(i)));
		}
	}
	
	WaveBasis::WaveBasis(double direction, freqtype k): direction_(direction), k_(k){}
	
	void WaveBasis::addBasis(freqtype k, pGeometry pgeom){
		pgeom->addBasis(pBasis(new WaveBasis(1.0, k)));
		pgeom->addBasis(pBasis(new WaveBasis(-1.0, k)));
	}

	WavePolBasis::WavePolBasis(int degree, double direction, freqtype k): degree_(degree), direction_(direction), k_(k){}
	
	void WavePolBasis::addBasis(int degree, freqtype k, pGeometry pgeom){
		for (int i=0;i<=degree;i++){
			pgeom->addBasis(pBasis(new WavePolBasis(i,1.0, k)));
			pgeom->addBasis(pBasis(new WavePolBasis(i,-1.0, k)));
		}
	}
	
	

}
