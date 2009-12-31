#ifndef QUADPOINTS_H_
#define QUADPOINTS_H_

#include<cstdlib>
#include "bem2d_defs.h"
#include "bem2d_exceptions.h"


namespace bem2d
{

struct QuadPoints1D {
        dvector x;
        dvector w;
        double a;
        double b;

};

struct QuadPoints2D {
        // Tensor Quadratur points in rectangle [a,b]x[c,d]
        // with coordinates (x[i],y[i]) and weights w[i]
        dvector x;
        dvector y;
        dvector w;
        double a;
        double b;
        double c;
        double d;

};


class Gauss1D
{
public:
        Gauss1D(int N) throw (ParameterException);
        // Create Gauss points between 0 and 1

        QuadPoints1D Rescale(double a, double b) const throw (IntervalException);

        inline const dvector& x() const {
                return x_;
        }

        inline const dvector& w() const {
                return w_;
        }

        inline const int Size() const {
                return N_;
        }

        ~Gauss1D();

private:
        dvector x_;
        dvector w_;
        int N_;
};

class Gauss2D
{
public:
        Gauss2D(int N) throw (ParameterException);
        // Create N*N Tensor-Gauss points in [0,1]^2

        QuadPoints2D Rescale(double a, double b, double c, double d) const throw (IntervalException);
        // Create Gauss Points in the rectangle [a,b]x[c,d]

        inline const dvector& x() const {
                return x_;
        }

        inline const dvector& y() const {
                return y_;
        }

        inline const dvector& w() const {
                return w_;
        }


        inline int Size() const {
                return N_; // Number of tensor Gauss points
        }

        ~Gauss2D();

private:
        dvector x_;
        dvector y_;
        dvector w_;
        int N_;
};

class AdaptedGauss1
{
        // Gauss Quadrature Points for a function that has singularities at the left side and lower side
        // of the box [0,1]x[0,1]
public:
        AdaptedGauss1(int N, int L, double sigma);

        inline const dvector& x() const {
                return x_;
        }

        inline const dvector& y() const {
                return y_;
        }

        inline const dvector& w() const {
                return w_;
        }


        inline int Size() const {
                return N_; // Number of tensor Gauss points
        }

        ~AdaptedGauss1();

private:
        dvector x_;
        dvector y_;
        dvector w_;
        int N_;
        int L_;
        double sigma_;
};

class AdaptedGauss2
{
        // Gauss Quadrature Points for a function that has singularities at the point (0,1)
public:
        AdaptedGauss2(int N, int L, double sigma);

        inline const dvector& x() const {
                return x_;
        }

        inline const dvector& y() const {
                return y_;
        }

        inline const dvector& w() const {
                return w_;
        }

        inline int Size() const {
                return N_; // Number of tensor Gauss points
        }

        ~AdaptedGauss2();

private:
        dvector x_;
        dvector y_;
        dvector w_;
        int N_;
        int L_;
        double sigma_;
};


class AdaptedGauss3
{
        // Gauss Quadrature Points for a function that has singularities at the point (1,0)
public:
        AdaptedGauss3(int N, int L, double sigma);

        inline const dvector& x() const {
                return x_;
        }

        inline const dvector& y() const {
                return y_;
        }

        inline const dvector& w() const {
                return w_;
        }

        inline int Size() const {
                return N_; // Number of tensor Gauss points
        }

        ~AdaptedGauss3();

private:
        dvector x_;
        dvector y_;
        dvector w_;
        int N_;
        int L_;
        double sigma_;
};


}
#endif /* QUADPOINTS_H_ */
