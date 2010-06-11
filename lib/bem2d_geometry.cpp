#include<algorithm>
#include "bem2d_geometry.h"


namespace bem2d
{


Geometry::Geometry(const std::vector<pElement>& elems):
                elements_(elems),
                size_(0)
{
        elements_bases_.resize(elems.size());

        // Now generate the point vector

        for (int i=0; i<elements_.size(); i++) {
                points_.push_back(elements_[i]->First());
        }

}
	Geometry::Geometry(std::vector<pGeometry>& geoms): size_(0)
{

    for (int i=0;i<geoms.size();i++){
        std::vector<pElement> elements=geoms[i]->elements();
        for (int j=0;j<elements.size();j++){
            elements_.push_back(elements[j]);
        }
    }
    int nelems=elements_.size();

    elements_bases_.resize(nelems);

    // Adjust the element ids and their connections
    
    for (int i=1;i<nelems-1;i++){
        elements_[i]->set_index(i);
        elements_[i]->set_next(i+1);
        elements_[i]->set_prev(i-1);
    }


    // Adjust first element

    elements_.front()->set_index(0);
    elements_.front()->set_next(1);
    elements_.front()->set_prev(nelems-1);

    // Adjust last element

    elements_.back()->set_index(nelems-1);
    elements_.back()->set_prev(nelems-2);
    elements_.back()->set_next(0);

    // Now generate the point vector

    for (int i=0; i<elements_.size(); i++) {
            points_.push_back(elements_[i]->First());
    }



}

}

