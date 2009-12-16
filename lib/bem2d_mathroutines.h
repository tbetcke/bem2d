#ifndef MATHROUTINES_H_
#define MATHROUTINES_H_

#include <vector>
#include<vector>
#include<cstdlib>
#include "bem2d_defs.h"
#include "bem2d_exceptions.h"
#include "bem2d_point.h"
#include "bem2d_geometry.h"


namespace bem2d
{
	
	
	// Special Functions
	
	complex BesselH0(double x);
	complex BesselH1(double x);
	
	// ----------------------------
	
	
	// Definition of matrix structures
	
	struct Matrix
	{
		explicit Matrix(int n);
		explicit Matrix(int m, int n);
		
		Matrix(const Matrix& m);
		pcvector data;
		int dim[2]; // M=dim[0]; N=dim[1];
		int msize; // Number of rows in local array
		int nsize; // Number of columns in local array	
#ifdef BEM2DMPI
		int desc[9];
#endif
	};
	
	typedef boost::shared_ptr<Matrix> pMatrix;
	
	
	Matrix operator+(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);
	Matrix operator-(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);
	
	Matrix operator*(const Matrix& lhs, complex& alpha);
	Matrix operator*(const complex& alpha, const Matrix& rhs);
	Matrix operator*(const Matrix& lhs, const double& alpha);
	Matrix operator*(const double& alpha, const Matrix& rhs);
	
	Matrix operator*(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);
	
	double L2Cond(const Matrix& m) throw (LapackError);
	
#ifdef BEM2DMPI	
	pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (ScaLapackError);
#else
	pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (LapackError);
#endif
	
	void InPolygon(const std::vector<pGeometry> polygons, const std::vector<Point>& testpoints,std::vector<int>& inpoly);
	void InPolygon(const pGeometry polygon, const std::vector<Point>& testpoints,std::vector<int>& inpoly);
	
	
}

#endif