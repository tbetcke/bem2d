#ifndef _SOUNDSOFTSCATTERING_H_
#define _SOUNDSOFTSCATTERING_H_

#include "bem2ddefs.h"
#include "geometry.h"
#include "Bem2dfun.h"
#include <vector>
#include "Outputhandler.h"
#include "quadrature.h"
#include "mathroutines.h"
#include "cblas.h"

namespace bem2d {
	
	template<typename T1, typename T2>
	class Soundsoftscattering {
	public:
		Soundsoftscattering(pGeometry geom, freqtype kvalue, T1 inc, T2 rhs);
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
		
		T1 incident;
		T2 frhs;
		
		pOutputhandler pout;
		
		QuadOption quadopts;
		matrix A;
		pcvector prhs;
		pcvector psol;
		
		
	};
	
	template<typename T1, typename T2>
	Soundsoftscattering<T1,T2>::Soundsoftscattering(pGeometry geom, freqtype kvalue, T1 inc, T2 rhs): pgeom(geom), 
	k(kvalue), A(geom->getsize()), frhs(rhs), incident(inc) {
		quadopts.L=3;
		quadopts.N=5;
		quadopts.sigma=0.15;
	}
	
	template<typename T1, typename T2>
	Soundsoftscattering<T1,T2>::~Soundsoftscattering(){}
	
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::setquadoption(int L, int N, double sigma){
		quadopts.L=L;
		quadopts.N=N;
		quadopts.sigma=sigma;
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::discretize(){
		
		
		bem2d::combinedsingleconjdouble scdl(k,frhs.geteta());
		A.data=discretekernel(*pgeom,quadopts,scdl);
		pcvector pident=dotprodbasfuns(*pgeom,quadopts,Idfun());
		cblas_zdscal(A.dim*A.dim,2.0,&(*A.data)[0],1);
		
		for (int i=0;i<A.dim;i++) (*A.data)[A.dim*i+i]+=(*pident)[i];
		
		prhs=bem2d::dotprodbasfuns(*pgeom,quadopts,frhs);
		cblas_zdscal(A.dim,2.0,&(*prhs)[0],1);
		
		/*
		std::cout << "A Matrix" << std::endl;
		for (int i=0;i<A.dim;i++) std::cout << (*A.data)[i] << std::endl;
		
		std::cout << "rhs";
		for (int i=0;i<A.dim;i++) std::cout << (*prhs)[i] << std::endl;
		std::cout << std::endl;
		*/
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::solve(){
		psol=pcvector(new cvector(*prhs));
		solve_system(A.data,psol);
		
		/*
		std::cout << "Sol" << std::endl;
		for (int i=0;i<A.dim;i++) std::cout << (*psol)[i] << std::endl;
		 */
	}
	
	template<typename T1,typename T2>
	pcvector Soundsoftscattering<T1,T2>::evalincident(const std::vector<Point>& points){
		pcvector pz(new cvector(points.size()));
#pragma omp for
		for (int i=0;i<points.size();i++) (*pz)[i]=incident(points[i]);
		return pz;
	}
	
	template<typename T1, typename T2>
	pcvector Soundsoftscattering<T1,T2>::evalsol(const std::vector<Point>& points){
		boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=pgeom->getflatmap();
		bem2d::singlelayer sl(k);
		pcvector kerneleval=pcvector(new cvector(points.size()));
		
		bem2d::Gauss1D g1d(quadopts.N);
		
		/*
		for (int i=0;i<A.dim;i++) {
			evalkernel((*bfuns)[i], points, *result, (*psol)[i], dl,g1d);
		}
		 */
		evalkernel(*bfuns,points,*kerneleval,*psol,sl,g1d);
		
		/*
		std::cout << "Single Layer Kernel";
		for (int i=0;i<points.size();i++) std::cout << (*kerneleval)[i] << std::endl;
		*/
		
		
		pcvector result=evalincident(points);
		complex alpha=-1.0;
		cblas_zaxpy(points.size(),&alpha, &(*kerneleval)[0], 1, &(*result)[0], 1);
		return result;
		
		 
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::setOutput(pOutputhandler output){
		pout=output;
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::writeIncident(){
		
		pcvector pinc=evalincident(pout->getmesh());
		pout->writeIncident(*pinc);
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::writeScattered(){
		pcvector psolution=evalsol(pout->getmesh());
		pcvector pincident=evalincident(pout->getmesh());
		pcvector pscattered=pcvector(new cvector(*psolution));
		complex alpha=-1.0;
		cblas_zaxpy(psolution->size(),&alpha, &(*pincident)[0], 1, &(*pscattered)[0], 1);
		pout->writeScattered(*pscattered);

	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::writeFull(){
		pcvector psolution=evalsol(pout->getmesh());
		pout->writeFull(*psolution);
	}
	
	template<typename T1, typename T2>
	void Soundsoftscattering<T1,T2>::writeAll(){
		pcvector pincident=evalincident(pout->getmesh());
		pcvector pfull=evalsol(pout->getmesh());
		
		pcvector pscattered=pcvector(new cvector(*pfull));
		complex alpha=-1.0;
		cblas_zaxpy(pincident->size(),&alpha, &(*pincident)[0], 1, &(*pscattered)[0], 1);
		pout->writeAll(*pincident,*pscattered,*pfull);
	}
	
	
	
	
}


#endif // _SOUNDSOFTSCATTERING_H_
