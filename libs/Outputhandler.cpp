#include "Outputhandler.h"
#include "outputroutines.h"

namespace bem2d{

	Gplotoutput::Gplotoutput(int xptsvalue, int yptsvalue, std::string namevalue): 
	xpts(xptsvalue),
	ypts(yptsvalue),
	name(namevalue){}
	
	void Gplotoutput::operator()(const std::vector<Point>& points, const cvector& vals){
		dvector realvals;
		realvals.reserve(vals.size());
		for (int i=0;i<vals.size();i++) realvals.push_back(vals[i].real());
		gplotout(name,points,realvals,xpts,ypts);
	}

}