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
	typedef double freqtype;
	
	
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
	}		
	
	
}

// External functions

extern "C"
{
	
	
	// Blacs Functions
	
	extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
	extern void   Cblacs_get( int context, int request, int* value);
	extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
	extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
	extern void   Cblacs_gridexit( int context);
	extern void   Cblacs_exit( int error_code);
	
	// Scalapack functions
	
	void sl_init_(int* ictxt, int* nprow, int* npcol);
	
	
}




#endif	/* _BEM2DEFS_H */

