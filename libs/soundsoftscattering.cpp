#include "soundsoftscattering.h"
#include "cblas.h"


namespace bem2d {


	Soundsoftscattering::Soundsoftscattering(pGeometry geom, freqtype kvalue): pgeom(geom), k(kvalue), A(geom->getsize()) {
		quadopts={3,3,0.15};
	}
	
	Soundsoftscattering::setincoming(pIncidentfun fun){
		incoming=fun;
	}
	
	Soundsoftscattering::setincoming(pIncidentfun fun, pIncidentfun normalderiv){
		incoming=fun;
		incomingnormal=normalderiv;
	}
	
	Soundsoftscattering::setquadoption(int L, int N, double sigma){
		quadopts.L=L;
		quadopts.N=N;
		quadopts.sigma=sigma;
	}
	
	void discretize(){
		
		bem2d::doublelayer g(k);
		A.data=bem2d::discretekernel(geom,opts,g);
		A=2.0*A+identity(A.dim);
		
		prhs=bem2d::discreterhs(geom,opts,pwave);
		cblas_zdscal(A.dim,2.0,&(*prhs)[0],1);
		
		
	}
		 
		
		
	};
	

}