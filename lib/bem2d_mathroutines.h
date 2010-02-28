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
complex BesselH0(complex x);
complex BesselH1(complex x);

// ----------------------------


// Definition of matrix structures

struct Matrix {
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

#ifdef BEM2DMPI
void L2NormCond(const Matrix& stiff, double& norm, double& cond) throw (ScaLapackError);
#else
void L2NormCond(const Matrix& stiff, double& norm, double& cond) throw (LapackError);
#endif

#ifdef BEM2DMPI
void SingularValues(const Matrix& stiff, pdvector& singvals) throw (ScaLapackError);
#else
void SingularValues(const Matrix& stiff, pdvector& singvals) throw (LapackError);
#endif


#ifdef BEM2DMPI
pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (ScaLapackError);
#else
pMatrix SolveSystem(Matrix& m, Matrix& rhs) throw (LapackError);
#endif

#ifdef BEM2DMPI
 void HermitianEigenvalues(const Matrix& k, pdvector& evalues) throw (ScaLapackError);
#else
 void HermitianEigenvalues(const Matrix& k, pdvector& evalues) throw (LapackError);
#endif

#ifdef BEM2DMPI
 void Eigenvalues(const Matrix& k, pcvector& evalues) throw (ScaLapackError);
#else
 void Eigenvalues(const Matrix& k, pcvector& evalues) throw (LapackError);
 void Eigenvalues(const Matrix& k, pcvector& evalues, pMatrix& evectors)
   throw (LapackError);
#endif

#ifdef BEM2DMPI
 Matrix ChangeBasis(const Matrix& k, const Matrix& m) throw (ScaLapackError);
#else
 Matrix ChangeBasis(const Matrix& k, const Matrix& m) throw (LapackError);
#endif

 void NumRange(const Matrix& m, int n, std::string filename);
 // Compute the boundary of the numerical range of a matrix m and write the points to a
 // file given by 'filename'. n is the number of discretization points from 0 to pi.

 Matrix ConjTranspose(const Matrix& m);

 Matrix ExtractColumn(const Matrix& m, int j);

 complex DotProduct(const Matrix& a, const Matrix& b) throw (ArrayMismatch);

void InPolygon(const std::vector<pGeometry> polygons, const std::vector<Point>& testpoints,std::vector<int>& inpoly);
void InPolygon(const pGeometry polygon, const std::vector<Point>& testpoints,std::vector<int>& inpoly);


}

#endif

