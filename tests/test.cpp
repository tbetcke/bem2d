#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>



class planewave {
public:
	planewave(bem2d::freqtype kvalue): k(kvalue){}
	bem2d::complex operator()(bem2d::Point x){
		return cos(k*x.x)+bem2d::complex(0,1)*sin(k*x.x);
	}
	void setnormal(bem2d::Point p){}
private:
	bem2d::freqtype k;
};
	
	

int main(int argc, char** argv){
	// File to test different things
	
	
	int n=500;
	
	bem2d::Trefoil tobj;
	bem2d::AnalyticCurve<bem2d::Trefoil> trefoil(n,tobj);
	
	bem2d::pGeometry pgeom=trefoil.GetGeometry();
	
	
	/*
	bem2d::Circle cobj;
	bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
	*/
	
	/*
	bem2d::DiskShapePiecewiseConst circle(n,1.0);
	*/
	/*
	bem2d::pGeometry pgeom=circle.GetGeometry();
	*/
	
	
	bem2d::freqtype k=20.0;
	bem2d::PolBasis::AddBasis(1,pgeom);
	std::cout << pgeom->size() << std::endl; 
	
	
	
	bem2d::PlaneWave pw(bem2d::Point(1,0),k);
	bem2d::CombinedPlaneWave cpw(bem2d::Point(1,0),k);
	
	
	std::cout << pw.k() << std::endl;

	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);

	soundsoft.Discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.Solve();

	std::cout << "System solved" << std::endl;
	
	
	int xpts=200; int ypts=200;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,2,-2,2,"trefoil"));
	soundsoft.SetOutput(pout);
	soundsoft.WriteAll();
	
	 
	 
	 
	return 0;
}
