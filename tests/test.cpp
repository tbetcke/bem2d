#include<iostream>
#include "../libs/bem2d.h"


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
	 */

	std::cout << std::endl;


	return 0;
}
