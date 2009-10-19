#include "geometry.h"

namespace bem2d {
	

	Geometry::Geometry(std::vector<pElement>& elem): 
	elements(elems),
	size(0) {
		basis_map.resize(elem.size());
	}
	
}

