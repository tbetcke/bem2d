#ifndef _BEM2DEFS_H
#define	_BEM2DEFS_H

#include<complex>
#include<vector>
#include "boost/shared_ptr.hpp"

// Some global typedefs, constants and external fct. definitions

namespace bem2d
{



// Type definitions

typedef std::complex<double> complex;
typedef std::vector<double> dvector;
typedef std::vector<complex> cvector;
struct freqtype
{
  double re;
  double im;
};


typedef boost::shared_ptr<cvector> pcvector;
typedef boost::shared_ptr<dvector> pdvector;




// Constants

const double PI=3.14159265358979323846;

extern "C" {
        void dstev_(char*, int*, double*, double*, double*, int*, double*, int*);
        void zaxpy_(int*, complex*, complex*, int*, complex*, int*);
        void zgetrf_(int*, int*, complex*, int*, int*, int*);
        void zgetrs_(char*, int*, int*, complex*, int*, int*, complex*, int*, int*);
        void zgesvd_(char*, char*, int*, int*, complex*, int*, double*, complex*, int*, complex*, int*, complex*, int*, double*, int*);
        void zpotrf_(const char*, const int*, complex*, const int*, int*);
  void zhegv_(const int* itype, const char* jobz, const char* uplo, const int* N, complex* A, 
	      const int* lda, complex* B, const int* ldb, double* w, complex* work, 
	      const int* lwork, double* rwork, int* info);
  void zheev_(const char* jobz, const char* uplo, const int* n, complex* A,
	      const int* lda, double* w, complex* work, const int* lwork,
	      double* rwork, int* info);
  void zgeev_(const char* jobvl, const char* jobvr, const int* n, 
	     complex* A, const int* lda, complex* w, complex* vl,
	     const int* ldvl, complex* vr, const int* ldvr, 
	     complex* work, const int* lwork, double* rwork,
	     int* info);
  void zbesh_(const double* zr, const double* zi,const double* fnu,
	      const int* kode, const int* M, const int* N, double* cyr,
	      double* cyi, int* nz, int* ierr);
}


}





#endif	/* _BEM2DEFS_H */

