#include<iostream>
#include "../libs/bem2d.h"
#include<cmath>


class outwave {
public:
	outwave(bem2d::freqtype kvalue): k(kvalue), s(kvalue){}
	bem2d::complex operator()(bem2d::Point x){
		bem2d::Point n(0,0);
		return s(n,x);
	}
	void setnormal(bem2d::Point p){}
private:
	bem2d::freqtype k;
	bem2d::singlelayer s;
};
	

int main(int argc, char** argv){
	// File to test different things

	
	/*
	std::cout << "Test" << std::endl;

	bem2d::dvector x;
	bem2d::dvector w;

	// bem2d::AdaptedGauss3 g(N,2,0.15);

	
	bem2d::Point p1(0,0);
	bem2d::Point p2(0,1);
	bem2d::Point p3(1,1);
	bem2d::Point p4(1,0);
	
	bem2d::pElement elem1(new bem2d::ConstElement(p1,p2,1));
	bem2d::pElement elem2(new bem2d::ConstElement(p2,p3,2));
	bem2d::pElement elem3(new bem2d::ConstElement(p3,p4,3));
	bem2d::pElement elem4(new bem2d::ConstElement(p4,p1,4));

	
	
	elem1->setNext(elem2->getIndex());
	elem2->setNext(elem3->getIndex());
	elem3->setNext(elem4->getIndex());
	elem4->setNext(elem1->getIndex());
	
	elem1->setPrev(elem4->getIndex());
	elem2->setPrev(elem1->getIndex());
	elem3->setPrev(elem2->getIndex());
	elem4->setPrev(elem3->getIndex());
	
	std::vector<bem2d::pElement> elements;
	elements.push_back(elem1);
	elements.push_back(elem2);
	elements.push_back(elem3);
	elements.push_back(elem4);
	
	bem2d::Geometry geom(elements);
	bem2d::pBasis b(new bem2d::ConstBasis);
	geom.addBasis(b);
	std::cout << geom.getsize() << std::endl; 
	
	boost::shared_ptr<std::vector<bem2d::complex> > pm;
	bem2d::singlelayer g(1.0);
	bem2d::QuadOption opts={5,3,0.15};
	pm=bem2d::discretekernel(geom,opts,g);
	
	
	int N=geom.getsize();
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			std::cout << (*pm)[N*j+i] << " ";
		}
		std::cout << std::endl;
	}
	
	boost::shared_ptr<std::vector<bem2d::complex> > prhs;
	outwave owave(1.0);
	prhs=bem2d::discreterhs(geom,opts,owave);
	
	std::cout << std::endl;
	
	for (int i=0;i<N;i++){
		std::cout << (*prhs)[i] << std::endl;
	}
	
	// Solve system
	
	bem2d::solve_system(pm,prhs);
	
	
	std::cout << std::endl;
	for (int i=0;i<N;i++){
		std::cout << (*prhs)[i] << std::endl;
	}
	
	
	
	
		/*
	double result=0;
	for (bem2d::dvector::const_iterator it=g.getpointsx();it!=g.pointendx();it++){
		std::cout <<  *(it)<< " ";
	}

	std::cout << std::endl;

	for (bem2d::dvector::const_iterator it=g.getpointsy();it!=g.pointendy();it++){
		std::cout <<  *(it)<< " ";
	}

	std::cout << std::endl;


	for (bem2d::dvector::const_iterator it=g.getweights();it!=g.weightend();it++){
		std::cout << *(it)<< " ";
		result +=*(it);
	}

	std::cout << std::endl;
	std::cout<< result <<std::endl;

	std::cout << std::endl;

	
	*/
	
	// Create a circular geometry
	
	
	int n=5; double k=1;
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
	
	bem2d::Geometry geom(elements);
	bem2d::pBasis b(new bem2d::ConstBasis);
	geom.addBasis(b);
	std::cout << geom.getsize() << std::endl; 
	
	boost::shared_ptr<std::vector<bem2d::complex> > pm;
	bem2d::doublelayer g(k);
	bem2d::QuadOption opts={5,3,0.15};
	pm=bem2d::discretekernel(geom,opts,g);
	
	
	boost::shared_ptr<std::vector<bem2d::complex> > prhs;
	outwave owave(k);
	prhs=bem2d::discreterhs(geom,opts,owave);

	// Update matrix and right-hand side
	
	int N=geom.getsize();
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			if (i==j){ 
				(*pm)[i*N+j]=1.0+2.0*(*pm)[i*N+j];
			}
			else {
				(*pm)[i*N+j]*=2.0;
			}
		}
		(*prhs)[i]*=2;
	}
	
	
	bem2d::solve_system(pm,prhs);
	
	// Now evaluate the solution on a grid
	
	int xpts=100; int ypts=100;
	bem2d::dvector x,y;
	boost::shared_ptr<std::vector<bem2d::Point> > pvec=bem2d::meshgrid(-2.0, 2.0,-2.0, 2.0, xpts, ypts);

	boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=geom.getflatmap();
	
	bem2d::cvector vals(xpts*ypts);
	bem2d::dvector realvals(xpts*ypts);
	bem2d::Gauss1D g1d(opts.N);
	
	
	for (int i=0;i<N;i++) {
		evalkernel((*bfuns)[i], *pvec, vals, (*prhs)[i], g,g1d);
		std::cout << vals[0] <<std::endl;
	}
	
	// Extract real part from vals array
	
	for (int i=0;i<N;i++) realvals[i]=vals[i].real();
	
	// Write solution to disk
	
	bem2d::gplotout("disk", *pvec, realvals, xpts, ypts);

	return 0;
}
