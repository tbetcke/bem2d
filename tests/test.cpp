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
	
	
	int n=1000;
	std::vector<bem2d::pElement> elements(n);
	std::vector<boost::shared_ptr<bem2d::Point> > points(n);
	
	double angle=2*M_PI/n;
	for (int i=0;i<n;i++) points[i]=boost::shared_ptr<bem2d::Point>(new bem2d::Point(cos(i*angle),sin(i*angle)));
	for (int i=0;i<n-1;i++) {
		elements[i]=bem2d::pElement(new bem2d::ConstElement(*(points[i]),*(points[i+1]),i));
	}
	elements[n-1]=bem2d::pElement(new bem2d::ConstElement(*(points[n-1]),*(points[0]),n-1));
	for (int i=1;i<n-1;i++) {
		elements[i]->setNext(i+1);
		elements[i]->setPrev(i-1);
	}
	elements[0]->setNext(1); elements[n-1]->setNext(0);
	
	bem2d::pGeometry pgeom(new bem2d::Geometry(elements));
	bem2d::pBasis b(new bem2d::ConstBasis);
	pgeom->addBasis(b);
	std::cout << pgeom->getsize() << std::endl; 
	
	
	
	bem2d::freqtype k=50;
	bem2d::PlaneWave pw(bem2d::Point(1,0),k);
	bem2d::Outwave owave(k);
	bem2d::Soundsoftscattering<bem2d::PlaneWave> soundsoft(pgeom,k,pw);
	//bem2d::pBem2dfun poutwave(new bem2d::Outwave(k));
	
	soundsoft.discretize();
	
	std::cout << "Solve system" << std::endl;

	soundsoft.solve();

	std::cout << "System solved" << std::endl;
	
	int xpts=400; int ypts=400;
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
