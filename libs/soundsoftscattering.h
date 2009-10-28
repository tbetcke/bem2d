#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include "bem2ddefs.h"
#include "Geometry.h"
#include "Bem2dfun.h"
#include <vector>
#include "Outputhandler.h"
#include "quadrature.h"
#include "mathroutines.h"
#include "cblas.h"

namespace bem2d {
	
	template<typename T>
	class Soundsoftscattering {
	public:
		Soundsoftscattering(pGeometry geom, freqtype kvalue, T fun);
		virtual ~Soundsoftscattering();
			
		
		void setquadoption(int L, int N, double sigma);
		virtual void discretize();
		
		void solve();
		
		pcvector evalincident(const std::vector<Point>& points);
		pcvector evalsol(const std::vector<Point>& points);
		
		void setOutput(pOutputhandler output);
		void writeIncident();
		void writeScattered();
		void writeFull();
		void writeAll();
		
	private:
		pGeometry pgeom;
		freqtype k;
		
		T incoming;
		
		pOutputhandler pout;
		
		QuadOption quadopts;
		matrix A;
		pcvector prhs;
		pcvector psol;
		
		
	};
	
	template<typename T>
	Soundsoftscattering<T>::Soundsoftscattering(pGeometry geom, freqtype kvalue, T fun): pgeom(geom), 
	k(kvalue), A(geom->getsize()), incoming(fun) {
		quadopts.L=3;
		quadopts.N=3;
		quadopts.sigma=0.15;
	}
	
	template<typename T>
	Soundsoftscattering<T>::~Soundsoftscattering(){}
	
	
	template<typename T>
	void Soundsoftscattering<T>::setquadoption(int L, int N, double sigma){
		quadopts.L=L;
		quadopts.N=N;
		quadopts.sigma=sigma;
	}
	
	template<typename T>
	void Soundsoftscattering<T>::discretize(){
		
		bem2d::doublelayer dl(k);
		A.data=discretekernel(*pgeom,quadopts,dl);
		pcvector pident=dotprodbasfuns(*pgeom,quadopts,Idfun());
		cblas_zdscal(A.dim*A.dim,2.0,&(*A.data)[0],1);
		
		for (int i=0;i<A.dim;i++) (*A.data)[A.dim*i+i]+=(*pident)[i];
		
		prhs=bem2d::dotprodbasfuns(*pgeom,quadopts,incoming);
		cblas_zdscal(A.dim,-2.0,&(*prhs)[0],1);
		
	}
	
	template<typename T>
	void Soundsoftscattering<T>::solve(){
		psol=pcvector(new cvector(*prhs));
		solve_system(A.data,psol);
	}
	
	template<typename T>
	pcvector Soundsoftscattering<T>::evalincident(const std::vector<Point>& points){
		pcvector pz(new cvector(points.size()));
#pragma omp for
		for (int i=0;i<points.size();i++) (*pz)[i]=incoming(points[i]);
		return pz;
	}
	
	template<typename T>
	pcvector Soundsoftscattering<T>::evalsol(const std::vector<Point>& points){
		boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=pgeom->getflatmap();
		bem2d::doublelayer dl(k);
		pcvector result=pcvector(new cvector(points.size()));
		
		bem2d::Gauss1D g1d(quadopts.N);
		
		/*
		for (int i=0;i<A.dim;i++) {
			evalkernel((*bfuns)[i], points, *result, (*psol)[i], dl,g1d);
		}
		 */
		evalkernel(*bfuns,points,*result,*psol,dl,g1d);
		return result;
	}
	
	template<typename T>
	void Soundsoftscattering<T>::setOutput(pOutputhandler output){
		pout=output;
	}
	
	template<typename T>
	void Soundsoftscattering<T>::writeIncident(){
		
		pcvector pinc=evalincident(pout->getmesh());
		pout->writeIncident(*pinc);
	}
	
	template<typename T>
	void Soundsoftscattering<T>::writeScattered(){
		pcvector psol=evalsol(pout->getmesh());
		pout->writeScattered(*psol);
	}
	
	template<typename T>
	void Soundsoftscattering<T>::writeFull(){
		pcvector pinc=evalincident(pout->getmesh());
		pcvector psol=evalsol(pout->getmesh());
		complex alpha(1.0,0);
		cblas_zaxpy(psol->size(),&alpha,&(*psol)[0],1,&(*pinc)[0],1);		
		pout->writeFull(*pinc);
	}
	
	template<typename T>
	void Soundsoftscattering<T>::writeAll(){
		pcvector pinc=evalincident(pout->getmesh());
		pcvector psol=evalsol(pout->getmesh());
		pcvector pfull(new cvector(pinc->size()));
		*pfull=*pinc;
		complex alpha=1.0;
		cblas_zaxpy(psol->size(),&alpha,&(*psol)[0],1,&(*pfull)[0],1);
		pout->writeAll(*pinc,*psol,*pfull);
	}
	
	
	
	
}


#endif // _SOUNDSOFTSCATTERING_H_
