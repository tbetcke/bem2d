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
	
	
	
	int n=1000; double k=50;
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
	planewave pwave(k);
	prhs=bem2d::discreterhs(geom,opts,pwave);

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
	
	int xpts=200; int ypts=200;
	bem2d::dvector x,y;
	boost::shared_ptr<std::vector<bem2d::Point> > pvec=bem2d::meshgrid(-2.0, 2.0,-2.0, 2.0, xpts, ypts);

	boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=geom.getflatmap();
	
	bem2d::cvector vals(xpts*ypts);
	bem2d::dvector realvals(xpts*ypts);
	bem2d::Gauss1D g1d(opts.N);
	
	
	for (int i=0;i<N;i++) {
		evalkernel((*bfuns)[i], *pvec, vals, (*prhs)[i], g,g1d);
	}
	
	// Extract real part from vals array
	
	for (int i=0;i<vals.size();i++) realvals[i]=vals[i].real();

	
	// Write solution to disk
	
	bem2d::gplotout("disk", *pvec, realvals, xpts, ypts);
	

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
