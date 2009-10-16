#ifndef _GEOMETRY_H
#define _GEOMETRY_H
#include<cstdlib>
#include<vector>
#include<map>


struct BasisComp {
	


class Geometry {
public:
	
private:
	std::vector<Element> elements;
	std::vector<Basis> bases;
	std::multimap<Element,Basis,BasisComp> elements_bases;
};


#endif // _GEOMETRY_H

