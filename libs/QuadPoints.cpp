#include<cmath>
#include<iostream>
#include "QuadPoints.h"
#include "quadrature.h"

namespace bem2d{

Gauss1D::Gauss1D(int N) throw (parameter_exception) {
	if (N<1) throw parameter_exception();
	this->N=N;
	gauss(x,w,N);
}

QuadPoints1D Gauss1D::rescale(double a, double b) const throw (interval_exception){
	if (a>=b) throw interval_exception();
	QuadPoints1D p;
	p.a=a; p.b=b; p.x=x; p.w=w;
	mappoints(p.x,p.w,a,b);
	return p;
}

Gauss1D::~Gauss1D() {
}

Gauss2D::Gauss2D(int N) throw (parameter_exception) {
	if (N<1) throw parameter_exception();

	int NN=N*N;
	this->N=NN;

	dvector xt; // 1d points
	dvector wt; // 1d weights
	gauss(xt, wt, N);

	x.resize(NN); y.resize(NN); w.resize(NN);

	for (int i=0; i<N;i++)
		for (int j=0;j<N;j++){
			x[i*N+j]=xt[i];
			y[i*N+j]=xt[j];
			w[i*N+j]=wt[i]*wt[j];
		}

}

QuadPoints2D Gauss2D::rescale(double a, double b, double c, double d) const throw (interval_exception){
	if (a>=b) throw interval_exception();
	if (c>=d) throw interval_exception();


	QuadPoints2D p;
	p.a=a; p.b=b; p.c=c; p.d=d;
	p.x=x; p.y=y; p.w=w;

	mappoints2d(p.x,p.y,p.w,a,b,c,d);
	return p;
}

Gauss2D::~Gauss2D() {
}

AdaptedGauss1::AdaptedGauss1(int N, int L, double sigma){

	this->N=N*N*L*L; this->L=L; this->sigma=sigma;

	dvector xl;
	for (int i=0;i<L;i++) xl.push_back(pow(sigma,L-i));
	xl.push_back(1.0);

	x.reserve(this->N); y.reserve(this->N); w.reserve(this->N);

	Gauss2D g2d(N);

	// Now iterate to create the right boxes and insert them


	for (int i=0;i<L;i++)
		for (int j=0;j<L;j++){
			QuadPoints2D q=g2d.rescale(xl[i],xl[i+1],xl[j],xl[j+1]);
			x.insert(x.end(),q.x.begin(),q.x.end());
			y.insert(y.end(),q.y.begin(),q.y.end());
			w.insert(w.end(),q.w.begin(),q.w.end());
		}

}

AdaptedGauss1::~AdaptedGauss1(){}


}
