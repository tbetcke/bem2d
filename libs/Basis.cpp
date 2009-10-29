#include "Basis.h"

namespace bem2d {
	
	Basis::~Basis(){}
	
	ConstBasis::ConstBasis(){
	};

	LegendrePolBasis::LegendrePolBasis(int degree): deg(degree){};

}
