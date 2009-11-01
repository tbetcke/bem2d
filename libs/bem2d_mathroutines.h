#ifndef MATHROUTINES_H_
#define MATHROUTINES_H_

#include <vector>
#include<vector>
#include<cstdlib>
#include "bem2d_defs.h"
#include "bem2d_exceptions.h"
#include "bem2d_point.h"


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
	
	matrix operator+(const matrix& lhs, const matrix& rhs) throw (ArrayMismatch);
	matrix operator-(const matrix& lhs, const matrix& rhs) throw (ArrayMismatch);
	
	matrix operator*(const matrix& lhs, complex& alpha);
	matrix operator*(const complex& alpha, const matrix& rhs);
	matrix operator*(const matrix& lhs, const double& alpha);
	matrix operator*(const double& alpha, const matrix& rhs);
			
	matrix operator*(const matrix& lhs, const matrix& rhs) throw (ArrayMismatch);
	
	void solve_system(pcvector pmatrix, pcvector prhs) throw (LapackError);
	
	
}

#endif