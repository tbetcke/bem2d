#ifndef MATHROUTINES_H_
#define MATHROUTINES_H_

#include<vector>
#include "bem2ddefs.h"
#include "exceptions.h"
#include<vector>
#include "Point.h"
#include<cstdlib>


namespace bem2d {
		
	
// Special Functions
	
	complex besselH0(double x);
	complex besselH1(double x);
	
// ----------------------------
	
	
// Definition of matrix structures
	
	struct matrix {
		matrix();
		explicit matrix(std::size_t n); 
		pcvector data;
		std::size_t dim;
	};
	
	struct identity: public matrix {
		identity(int n);
	};
	
	matrix operator+(const matrix& lhs, const matrix& rhs) throw (array_mismatch);
	matrix operator-(const matrix& lhs, const matrix& rhs) throw (array_mismatch);
	
	matrix operator*(const matrix& lhs, complex& alpha);
	matrix operator*(const complex& alpha, const matrix& rhs);
	matrix operator*(const matrix& lhs, const double& alpha);
	matrix operator*(const double& alpha, const matrix& rhs);
			
	matrix operator*(const matrix& lhs, const matrix& rhs) throw (array_mismatch);
	
	void solve_system(pcvector pmatrix, pcvector prhs) throw (lapack_error);
	
	
}

#endif