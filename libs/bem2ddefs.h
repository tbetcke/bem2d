#ifndef _BEM2DEFS_H
#define	_BEM2DEFS_H

#include<complex>
#include<vector>

// Typedefs and definitions of external functions

namespace bem2d {


    // Type definitions

    typedef std::complex<double> complex;
    typedef std::vector<double> dvector;
	typedef double freqtype;

    // External functions

    extern "C" {

        void dstev_(char*, int*, double*, double*, double*, int*, double*, int*);
        void zaxpy_(int*, complex*, complex*, int*, complex*, int*);
        void zgetrf_(int*, int*, complex*, int*, int*, int*);
        void zgetrs_(char*, int*, int*, complex*, int*, int*, complex*, int*, int*);
    }



}



#endif	/* _BEM2DEFS_H */

