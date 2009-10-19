#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include<cstdlib>
#include<vector>
#include<map>
#include<boost/shared_ptr.hpp>

namespace bem2d{

	
class Geometry {
public:	

	typedef std::vector< std::map<std::size_t,pBasis> > basis_map;
	typedef std::vector< std::pair<pElement,pBasis> > flat_basis_map;
	
	Geometry(const std::vector<pElement>& elems);
	

	inline const std::vector<pElement>& getElements(){
		return elements;
	}
	
	inline void addBasis(pBasis b){
		for (std::size_t i=0; i<elements.size();i++){
			elements_bases[i].insert(std::pair<std::size_t,pBasis>(size,b));
			size++;
		}
	}
		
	inline void addBasis(std::size_t i,pBasis b){
		elements_bases[i].insert(std::pair<std::size_t,pBasis>(size,b));
		size++;
	}
	
	inline boost::shared_ptr<flat_basis_map> getflatmap() const {
		boost::shared_ptr<flat_basis_map> p(new flat_basis_map);
		for (int i=0;i<elements_bases.size();i++){
			std::map<std::size_t,pBasis>::const_iterator jt;
			for (jt=elements_bases[i].begin(); jt!=elements_bases[i].end();jt++){
				p->push_back(std::pair<pElement,pBasis>(elements[i],jt->second));
			}
		}
		return p;
	}
		
	
private:
	std::vector<pElement> elements;
	basis_map elements_bases;
	int size; // Number of basis functions
};

}

#endif // _GEOMETRY_H

