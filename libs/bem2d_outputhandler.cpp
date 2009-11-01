#include "boost/shared_ptr.hpp"
#include "bem2d_outputhandler.h"
#include "bem2d_outputroutines.h"

namespace bem2d{

	pdvector Outputhandler::realdata(const cvector& vals){
		pdvector result(new dvector(vals.size()));
		for (int i=0;i<vals.size();i++) (*result)[i]=vals[i].real();
		return result;
	}
	
	pdvector Outputhandler::imagdata(const cvector& vals){
		pdvector result(new dvector(vals.size()));
		for (int i=0;i<vals.size();i++) (*result)[i]=vals[i].imag();
		return result;
	}
	
	pdvector Outputhandler::turntorealimag(const cvector& vals){
		if (real){
			return realdata(vals);
		}
		else {
			return imagdata(vals);
		}
	}
	

	Gplotoutput::Gplotoutput(int xptsvalue, int yptsvalue, int ax, int bx, int ay, int by, std::string fname){
		xpts=xptsvalue;
		ypts=yptsvalue;
		name=fname;
		setreal();
		mesh=*(meshgrid(ax,bx,ay,by,xpts,ypts));
	}
	
	void Gplotoutput::writeIncident(const cvector& vals){
		gplotout(name+"_inc",mesh,*(turntorealimag(vals)),xpts,ypts);
	}
	
	void Gplotoutput::writeScattered(const cvector& vals){
		gplotout(name+"scatt",mesh,*(turntorealimag(vals)),xpts,ypts);
	}
	
	void Gplotoutput::writeFull(const cvector& vals){
		gplotout(name+"full",mesh,*(turntorealimag(vals)),xpts,ypts);		
	}
	
	void Gplotoutput::writeAll(const cvector& valsincident, const cvector& valsscattered, const cvector& valsfull){
		gplotout(name+"_inc",mesh,*(turntorealimag(valsincident)),xpts,ypts);
		gplotout(name+"scatt",mesh,*(turntorealimag(valsscattered)),xpts,ypts);
		gplotout(name+"full",mesh,*(turntorealimag(valsfull)),xpts,ypts);		
	}
			
	
}

