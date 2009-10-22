#include "geometry.h"

namespace bem2d {
	

	Geometry::Geometry(const std::vector<pElement>& elems): 
	elements(elems),
	size(0) {
		elements_bases.resize(elems.size());
	}
	
}

