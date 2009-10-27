#include "mathroutines.h"
#include "gsl/gsl_sf_bessel.h"
#include "fstream"
#include "cblas.h"

namespace bem2d {
	

	complex besselH0(double x) {
        double fr, fi;
		
        // Evaluate Bessel function
		
        fr = gsl_sf_bessel_J0(x);
        fi = gsl_sf_bessel_Y0(x);
        return complex(fr, fi);
		
    }
	
    complex besselH1(double x) {
        double fr, fi;
		
        // Evaluate Bessel function
		
        fr = gsl_sf_bessel_J1(x);
        fi = gsl_sf_bessel_Y1(x);
        return complex(fr, fi);
		
    }
	
	
	matrix::matrix(): dim(0){
		data=pcvector(new cvector);
	
	};
	
	matrix::matrix(std::size_t n): dim(n){
		data=pcvector(new cvector);
		data->resize(n*n);
	}
	
	identity::identity(int n): matrix(n){
		for (int i=0;i<n;i++) (*data)[n*i+i]=1.0;
	}
	
	
	matrix operator+(const matrix& lhs, const matrix& rhs) throw (array_mismatch){
		if (lhs.dim!=rhs.dim) throw array_mismatch();
		matrix result(lhs.dim);
		*result.data=*rhs.data;
		complex alpha=1.0;
		cblas_zaxpy(lhs.dim*lhs.dim,&alpha,&(*lhs.data)[0],1,&(*result.data)[0],1);
		return result;
		
	}
	
	matrix operator-(const matrix& lhs, const matrix& rhs) throw (array_mismatch){
		if (lhs.dim!=rhs.dim) throw array_mismatch();
		matrix result(lhs.dim);
		*result.data=*lhs.data;
		complex alpha=-1.0;
		cblas_zaxpy(lhs.dim*lhs.dim,&alpha,&(*rhs.data)[0],1,&(*result.data)[0],1);
		return result;
		
	}
	
	
	matrix operator*(const matrix& lhs, const complex& alpha){
		matrix result(lhs.dim);
		*result.data=*lhs.data;
		cblas_zscal(lhs.dim*lhs.dim,&alpha,&(*result.data)[0],1);
		return result;
	}
	
	matrix operator*(const complex& alpha, const matrix& rhs){
		return rhs*alpha;
	}
	
	matrix operator*(const matrix& lhs, const double& alpha){
		matrix result(lhs.dim);
		*result.data=*lhs.data;
		cblas_zdscal(lhs.dim*lhs.dim,alpha,&(*result.data)[0],1);
		return result;
	}
				
	matrix operator*(const double& alpha, const matrix& rhs){
		return rhs*alpha;
	}
	
	matrix operator*(const matrix& lhs, const matrix& rhs) throw (array_mismatch){
		if (lhs.dim!=rhs.dim) throw array_mismatch();
		matrix result(lhs.dim);
		complex alpha=1.0;
		complex beta=0.0;
		cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,lhs.dim,rhs.dim,lhs.dim,&alpha,&(*lhs.data)[0],lhs.dim,&(*rhs.data)[0],
					rhs.dim,&beta,&(*result.data)[0],rhs.dim);
		return result;
		
	}
	
	
	
	void solve_system(pcvector pmatrix, pcvector prhs) throw (lapack_error){
		
		int N=prhs->size();
		int info;
		char trans ='N';
		int nrhs=1;
		std::vector<int> ipiv(N);
		
		
		zgetrf_(&N, &N, &((*pmatrix)[0]), &N, &(ipiv[0]),&info);
		if (info!=0) throw lapack_error();
		
		zgetrs_(&trans,&N,&nrhs,&((*pmatrix)[0]),&N,&(ipiv[0]),&((*prhs)[0]),&N,&info);
		if (info!=0) throw lapack_error();
		
	}
	
	
	
	
	
}
