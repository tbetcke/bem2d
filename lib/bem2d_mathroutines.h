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
  Matrix();
  explicit Matrix(std::size_t n);
  Matrix(const Matrix& m);
  pcvector data;
  std::size_t dim;
};

struct Identity: public Matrix
{
  Identity(int n);
};

Matrix operator+(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);
Matrix operator-(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);

Matrix operator*(const Matrix& lhs, complex& alpha);
Matrix operator*(const complex& alpha, const Matrix& rhs);
Matrix operator*(const Matrix& lhs, const double& alpha);
Matrix operator*(const double& alpha, const Matrix& rhs);

Matrix operator*(const Matrix& lhs, const Matrix& rhs) throw (ArrayMismatch);
	
	double L2Cond(const Matrix& m) throw (LapackError);

void SolveSystem(pcvector pmatrix, pcvector prhs) throw (LapackError);

void InPolygon(const std::vector<pGeometry> polygons, const std::vector<Point>& testpoints,std::vector<int>& inpoly);
void InPolygon(const pGeometry polygon, const std::vector<Point>& testpoints,std::vector<int>& inpoly);

	
}

#endif