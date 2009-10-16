#include<iostream>
#include "../libs/bem2d.h"


int main(int argc, char** argv){
	// File to test different things

	std::cout << "Test" << std::endl;

	bem2d::dvector x;
	bem2d::dvector w;
	int N=3;

	bem2d::AdaptedGauss3 g(N,2,0.15);


	
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


	return 0;
}
