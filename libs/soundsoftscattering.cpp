#include "soundsoftscattering.h"
#include "kernel.h"
#include "cblas.h"


namespace bem2d {


	Soundsoftscattering::Soundsoftscattering(pGeometry geom, freqtype kvalue): pgeom(geom), k(kvalue), A(geom->getsize()) {
		quadopts.L=3;
		quadopts.N=3;
		quadopts.sigma=0.15;
	}
	
	Soundsoftscattering::~Soundsoftscattering(){}
	
	void Soundsoftscattering::setincoming(pBem2dfun fun){
		incoming=fun;
	}
	
	void Soundsoftscattering::setincoming(pBem2dfun fun, pBem2dfun normalderiv){
		incoming=fun;
		incomingnormal=normalderiv;
	}
	
	void Soundsoftscattering::setquadoption(int L, int N, double sigma){
		quadopts.L=L;
		quadopts.N=N;
		quadopts.sigma=sigma;
	}
	
	void Soundsoftscattering::discretize(){
		
		bem2d::doublelayer dl(k);
		A.data=discretekernel(*pgeom,quadopts,dl);
		pBem2dfun id(new Idfun());
		pcvector pident=dotprodbasfuns(*pgeom,quadopts,*id);
		
		for (int i=0;i<A.dim;i++) (*A.data)[A.dim*i+i]+=(*pident)[i];
		
		
		prhs=bem2d::dotprodbasfuns(*pgeom,quadopts,*incoming);
		cblas_zdscal(A.dim,-2.0,&(*prhs)[0],1);
			
		
	}
	
	void Soundsoftscattering::solve(){
		psol=pcvector(new cvector(*prhs));
		solve_system(A.data,psol);
	}
	
	pcvector Soundsoftscattering::evalincident(const std::vector<Point>& points){
		pcvector pz(new cvector(points.size()));
		for (std::size_t i=0;i<points.size();i++) (*pz)[i]=(*incoming)(points[i]);
		return pz;
	}
		
	
	pcvector Soundsoftscattering::evalsol(const std::vector<Point>& points){
		boost::shared_ptr<bem2d::Geometry::flat_basis_map> bfuns=pgeom->getflatmap();
		bem2d::doublelayer dl(k);
		pcvector result=pcvector(new cvector(points.size()));
		
		bem2d::Gauss1D g1d(quadopts.N);			
		for (int i=0;i<A.dim;i++) {
			evalkernel((*bfuns)[i], points, *result, (*psol)[i], dl,g1d);
		}
		return result;
	}
		
	void Soundsoftscattering::setOutput(pOutputhandler output){
		pout=output;
	}
	
	void Soundsoftscattering::plotIncident(const std::vector<Point>& points){
		pcvector pinc=evalincident(points);
		(*pout)(points,*pinc);
	}
	
	void Soundsoftscattering::plotScattered(const std::vector<Point>& points){
		pcvector psol=evalsol(points);
		(*pout)(points,*psol);
	}
	
	void Soundsoftscattering::plotFull(const std::vector<Point>& points){
		pcvector psol=evalsol(points);
		pcvector pinc=evalincident(points);
		complex alpha(1.0,0);
		cblas_zaxpy(points.size(),&alpha,&(*psol)[0],1,&(*pinc)[0],1);		
		(*pout)(points,*pinc);
	}		
		
		
}