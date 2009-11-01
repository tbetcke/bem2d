#include<iostream>
#include "../libs/bem2d.h"
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
	
	//bem2d::diskshape_piecewise_const disk(n,1.0);
		
	/*
	bem2d::trefoil tobj;
	bem2d::analytic_curve<bem2d::trefoil> trefoil(n,tobj);
	
	bem2d::pGeometry pgeom=trefoil.getGeometry();
	*/
	
	
	bem2d::Circle cobj;
	bem2d::analytic_curve<bem2d::Circle> circle(n,cobj);
	
	bem2d::pGeometry pgeom=circle.getGeometry();
	
	
	/* 
	bem2d::pBasis b0(new bem2d::PolBasis(0));
	bem2d::pBasis b1(new bem2d::PolBasis(1));
	bem2d::pBasis b2(new bem2d::PolBasis(2));
	pgeom->addBasis(b0);
	pgeom->addBasis(b1);
	pgeom->addBasis(b2);
	 */
	bem2d::freqtype k=10;
	bem2d::PolBasis::addBasis(0,pgeom);
	std::cout << pgeom->size() << std::endl; 
	
	
	
	bem2d::PlaneWave pw(bem2d::Point(1,0),k);
	bem2d::CombinedPlaneWave cpw(bem2d::Point(1,0),k);

	//bem2d::Outwave owave(k);
	
	bem2d::Soundsoftscattering<bem2d::PlaneWave,bem2d::CombinedPlaneWave> soundsoft(pgeom,k,pw,cpw);
	
	soundsoft.discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.solve();

	std::cout << "System solved" << std::endl;
	
	
	int xpts=100; int ypts=100;
	bem2d::pOutputhandler pout(new bem2d::Gplotoutput(xpts,ypts,-2,2,-2,2,"disk"));
	soundsoft.setOutput(pout);
	soundsoft.writeAll();
	
	 
	 
	/*
	bem2d::singlelayer sl(1.0);
	std::cout << sl(bem2d::Point(0,0),bem2d::Point(-2,-2)) << std::endl;
	
	std::vector<bem2d::Point> pz(1); pz[0]=bem2d::Point(-2,-2);
	bem2d::pcvector fvals=soundsoft.evalsol(pz);
	
	std::cout << (*fvals)[0] << std::endl;
	*/
	
	// for (int i=0;i<psol->size();i++) std::cout << (*psol)[i] << std::endl;
	// std::cout << std::endl;
	
	

	/*
	bem2d::identity a(2);
	bem2d::matrix b(2);
	(*b.data)[2]=5;
	(*b.data)[0]=6;
	
	bem2d::matrix c=a*b;
	
	for (int i=0;i<c.dim;i++){
		for (int j=0;j<c.dim;j++){
			std::cout << (*c.data)[i+c.dim*j] << " ";
		}
		std::cout << std::endl;
	}
	
	std::cout << std::endl;
	*/
	 
	return 0;
}
