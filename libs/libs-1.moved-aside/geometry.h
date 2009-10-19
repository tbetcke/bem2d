#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include<cstdlib>
#include<vector>
#include<map>

namespace bem2d{

struct ElementComp {
	
	bool operator()(const pElement& lhs, const pElement& rhs) const {
		return lhs->getIndex()<rhs->getIndex();
	}
};
	
class Geometry {
public:	

	typedef std::vector< std::map<std::size_t,pElement> > basis_map;
	
	inline void addElement(pElement elem){
		elements.insert(elements.end(),elem);
	}
	
	inline void addElement(const std::vector<pElement>& elems){
		elements.insert(elements.end(),elems.begin(),elems.end());
	}


	inline const std::vector<pElement>& getElements(){
		return elements;
	}
	
	inline const std::multimap<pElement,pBasis,ElementComp>& getElemMap(){
		return elements_bases;
	}
	
	inline void addBasis(pBasis b){
		for (std::vector<pElement>::const_iterator it=elements.begin(); it!=elements.end();it++){
			elements_bases.insert(std::pair<pElement,pBasis>(*it,b));
		}
	}
		
	inline void addBasis(pElement elem,pBasis b){
		elements_bases.insert(std::pair<pElement,pBasis>(elem,b));
	}
	
	
private:
	std::vector<pElement> elements;
	std::multimap<pElement,pBasis,ElementComp> elements_bases;
};

}

#endif // _GEOMETRY_H

