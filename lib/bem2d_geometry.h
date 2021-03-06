#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include <cstdlib>
#include <vector>
#include <map>
#include "boost/shared_ptr.hpp"
#include "bem2d_basis.h"
#include "bem2d_element.h"
#include "bem2d_point.h"

namespace bem2d
{

class Geometry;
typedef boost::shared_ptr<Geometry> pGeometry;

class Geometry
{
public:

        typedef std::vector< std::map<std::size_t,pBasis> > basis_map;
        typedef std::vector< std::pair<pElement,pBasis> > flat_basis_map;

        Geometry(const std::vector<pElement>& elems);
        Geometry(std::vector<pGeometry>& geoms);

        inline std::vector<pElement>& elements() {
                return elements_;
        }

	inline const Geometry::basis_map elements_bases(){
	  return elements_bases_;
	}

        inline void AddBasis(pBasis b) {
		std::cout << "Elems: " << elements_.size() << std::endl;
                for (std::size_t i=0; i<elements_.size(); i++) {
                        elements_bases_[i].insert(std::pair<std::size_t,pBasis>(size_,b));
                        size_++;
                }
        }

        inline void AddBasis(std::size_t i,pBasis b) {
                elements_bases_[i].insert(std::pair<std::size_t,pBasis>(size_,b));
                size_++;
        }

        inline boost::shared_ptr<flat_basis_map> FlatMap() const {
                boost::shared_ptr<flat_basis_map> p(new flat_basis_map);
                for (int i=0; i<elements_bases_.size(); i++) {
                        std::map<std::size_t,pBasis>::const_iterator jt;
                        for (jt=elements_bases_[i].begin(); jt!=elements_bases_[i].end(); jt++) {
                                p->push_back(std::pair<pElement,pBasis>(elements_[i],jt->second));
                        }
                }
                return p;
        }


        inline std::size_t size() const {
                return size_;
        }


        inline const PointVector& points() const {
                return points_;
        }


private:
        std::vector<pElement> elements_;
        basis_map elements_bases_;
        PointVector points_; // Stores the points defining a geometry
        int size_; // Number of basis functions
};


}

#endif // _GEOMETRY_H

