#include<iostream>
#include "../lib/bem2d.h"
#include<cmath>
#include<ctime>


int main(int argc, char** argv){
	
    int m=2;
	bem2d::freqtype k=5*m;
	int n=20*m;
	
	double a=0.31;
	double c=1;
	double l=c-a;
	
	
	
	std::vector<bem2d::Point> trapping;
	trapping.push_back(bem2d::Point(0,0));
	trapping.push_back(bem2d::Point(-c,0));
	trapping.push_back(bem2d::Point(-c,-l));
	trapping.push_back(bem2d::Point(l,-l));
	trapping.push_back(bem2d::Point(l,2*c-l));
	trapping.push_back(bem2d::Point(-c,2*c-l));
	trapping.push_back(bem2d::Point(-c,2*a));
	trapping.push_back(bem2d::Point(0,2*a));
	bem2d::Polygon poly(trapping,n);
	bem2d::pGeometry pgeom=poly.GetGeometry();
	
	
	
	bem2d::PolBasis::AddBasis(0,pgeom);
	
	
	bem2d::Point direction=bem2d::normalize(bem2d::Point(1,0));
	bem2d::PlaneWave pw(direction,k);
	bem2d::CombinedPlaneWave cpw(direction,k);
	
	bem2d::SoundSoftScattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	soundsoft.SetQuadOption(5,5,0.15);
	soundsoft.set_polygons(pgeom);
	soundsoft.set_plotInterior();

	
	
	std::cout << "Discretizing with k=" << k << " and n=" << pgeom->size() << std::endl; 
	clock_t start, finish;
	double time;
	
	start=clock();
	soundsoft.Discretize();
	finish=clock();
	time=(double(finish)-double(start))/CLOCKS_PER_SEC;
	std::cout << "Computing time (minutes): " << time/60 << std::endl;
	std::cout << "Condition Number: " << soundsoft.L2Condition() << std::endl;
	
	int xpts=100; int ypts=100;
	bem2d::pOutputHandler pout(new bem2d::GplotOutput(xpts,ypts,-2,3,-2,2,"trapping"));
	soundsoft.SetOutput(pout);	
	soundsoft.WriteAll();

	
	
	
	
	
	return 0;
}



