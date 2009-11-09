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
	
	
	int n=400;
	
	std::vector<bem2d::Point> square;
	square.push_back(bem2d::Point(-.5,-.5));
	square.push_back(bem2d::Point(.5,-.5));
	square.push_back(bem2d::Point(.5,.5));
	square.push_back(bem2d::Point(-.5,.5));
	bem2d::Polygon poly(square,n);
	bem2d::pGeometry pgeom=poly.GetGeometry();
	
	/*
	bem2d::Trefoil tobj;
	bem2d::AnalyticCurve<bem2d::Trefoil> trefoil(n,tobj);
	
	bem2d::pGeometry pgeom=trefoil.GetGeometry();
	*/
	
	/*
	bem2d::Circle cobj;
	bem2d::AnalyticCurve<bem2d::Circle> circle(n,cobj);
	bem2d::pGeometry pgeom=circle.GetGeometry();
	*/
	
	
	/*
	bem2d::DiskShapePiecewiseConst circle(n,1.0);
	*/
	
	
	bem2d::freqtype k=20;
	bem2d::PolBasis::AddBasis(0,pgeom);
	std::cout << pgeom->size() << std::endl; 
	
	
	bem2d::Point direction=bem2d::normalize(bem2d::Point(1,-1));
	bem2d::PlaneWave pw(direction,k);
	bem2d::CombinedPlaneWave cpw(direction,k);
	
	
	std::cout << pw.k() << std::endl;

	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);
	soundsoft.set_polygons(pgeom);
	soundsoft.set_plotInterior();

	soundsoft.Discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.Solve();

	std::cout << "System solved" << std::endl;
	
	std::cout << "Condition Number: " << soundsoft.L2Condition() << std::endl;
	
	int xpts=100; int ypts=100;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-1,1,-1,1,"square"));
	soundsoft.SetOutput(pout);
	soundsoft.WriteAll();
	
	
	
	 
	 
	 
	return 0;
}
