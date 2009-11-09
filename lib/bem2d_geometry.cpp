#include "bem2d_geometry.h"

namespace bem2d
{


Geometry::Geometry(const std::vector<pElement>& elems):
    elements_(elems),
    size_(0)
{
  elements_bases_.resize(elems.size());
	
  // Now generate the point vector
	
	for (int i=0;i<elements_.size(); i++){
		points_.push_back(elements_[i]->First());
	}
	
}

}

