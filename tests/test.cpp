#include<iostream>
#include "../libs/bem2d.h"


int main(int argc, char** argv){
	// File to test different things

	std::cout << "Test" << std::endl;

	bem2d::dvector x;
	bem2d::dvector w;
	int N=50;

	bem2d::gauss(x,w,N);


	for (int i=0;i<N;i++) std::cout << x[i]<< std::endl;
	std::cout << std::endl;
	for (int i=0;i<N;i++) std::cout << w[i]<< std::endl;

	return 0;
}
